/*
 * $Id$ ::2B6A7A!
 * Copyright (c) 2005 Dr. Vahan Simonyan and Dr. Raja Mazumder.
 * This software is protected by U.S. Copyright Law and International
 * Treaties. Unauthorized use, duplication, reverse engineering, any
 * form of redistribution, use in part or as a whole, other than by
 * prior, express, written and signed agreement is subject to penalties.
 * If you have received this file in error, please notify copyright
 * holder and destroy this and any other copies. All rights reserved.
 */
#include <qlib/QPrideProc.hpp>
#include <slib/utils.hpp>
#include <ssci/bio.hpp>
#include <violin/violin.hpp>
#include <qpsvc/qpsvc-dna-alignment-remapper.hpp>

#include "bioclones.hpp"
#include "cloneTrees.hpp"

#define _getPairedRead(_iqry) (_iqry + (1-(2*(_iqry/(_bioal_long_cnt/2))))*(_bioal_long_cnt/2))
#define REALFNLNM "re-alignment.hiveal"
#define POPULSVCNAME "dna-popul"

class DnaPopul: public sQPrideProc {
public:
	DnaPopul(const char * defline00, const char * srv) :
			sQPrideProc(defline00, srv) {
	}

	idx launchRemappers(sHiveal & hiveal);
	idx constructPairMap(sVec<idx> &alignPairMap, sVec<idx> &alignInd,
			sBioal * bioal, const char * query);

	virtual idx OnExecute(idx);
};

idx DnaPopul::constructPairMap(sVec<idx> &alignPairMap, sVec<idx> &alignInd,
		sBioal * bioal, const char * query) {

	idx _bioal_long_cnt = bioal->Qry->getlongCount();
	if (_bioal_long_cnt % 2)
		return 0;

	sVec<idx> reads_short2long(sMex::fSetZero | sMex::fExactSize);
	reads_short2long.resize(_bioal_long_cnt);

	idx shortI = 0;
	sVec<idx> reads_ind_short2long(sMex::fSetZero | sMex::fExactSize);
	sVec<idx> reads_cnt_short2long(sMex::fSetZero | sMex::fExactSize);
	reads_ind_short2long.resize(bioal->Qry->dim());
	reads_cnt_short2long.resize(bioal->Qry->dim());
	idx cnt = 0;
	for (idx i = 0; i < bioal->Qry->dim(); ++i) {
		reads_ind_short2long[i] = cnt;
		cnt += bioal->Qry->rpt(i);
	}

	idx curInd = 0;
	progress100Start = 1;
	progress100End = 3;
	for (idx i = 0; i < _bioal_long_cnt; ++i) {

		if (!reqProgress(1, i, _bioal_long_cnt))
			return 0;

		shortI = bioal->Qry->long2short(i);
		curInd = reads_ind_short2long[shortI] + reads_cnt_short2long[shortI]++;
		reads_short2long[curInd] = i;
	}
	sDic<idx> long2al;
	alignPairMap.cut(0);
	alignPairMap.setflag(sMex::fSetZero | sMex::fExactSize);
	alignPairMap.resize(_bioal_long_cnt);
	alignInd.cut(0);
	alignInd.setflag(sMex::fSetZero | sMex::fExactSize);
	alignInd.resize(2 * bioal->dimAl());
	idx qry = 0, iLqry = 0, rptcnt = 0;
	cnt = 0;
	progress100Start = 3;
	progress100End = 6;
	for (idx ia = 0; ia < bioal->dimAl(); ++ia) {
		if (!reqProgress(1, ia, bioal->dimAl()))
			return 0;

		qry = bioal->getAl(ia)->idQry();
		rptcnt = (
				(qry + 1 >= reads_ind_short2long.dim()) ?
						reads_short2long.dim() : reads_ind_short2long[qry + 1])
				- reads_ind_short2long[qry];

		for (idx iq = 0; iq < rptcnt; ++iq) {
			iLqry = reads_short2long[reads_ind_short2long[qry] + iq];
			*long2al.set(&iLqry, sizeof(iLqry)) = ia;
		}
	}
	idx * paired_ial = 0, ind_cnt = 0, pair_al_cnt = 0, ia_cnt = 0, paired_qry =
			0, long_qry;
	progress100Start = 6;
	progress100End = 9;

	for (idx ia = 0; ia < bioal->dimAl(); ++ia) {

		if (!reqProgress(1, ia, bioal->dimAl()))
			return 0;

		qry = bioal->getAl(ia)->idQry();
		ia_cnt = reads_cnt_short2long[qry];
		for (idx iq = 0; iq < ia_cnt; ++iq) {
			long_qry = reads_short2long[reads_ind_short2long[qry] + iq];
			paired_qry = _getPairedRead(long_qry);
			paired_ial = long2al.get(&paired_qry, sizeof(paired_qry));
			if (paired_ial) {
				if (!alignInd[2 * ia]) {
					alignInd[2 * ia + 1] = ind_cnt;
					++pair_al_cnt;
				}
				alignInd[2 * ia]++;
				alignPairMap[ind_cnt++] = *paired_ial;
			}
		}
	}
	progress100Start = 10;
	progress100End = 100;
	alignPairMap.cut(ind_cnt);
	return pair_al_cnt;
}

idx DnaPopul::launchRemappers(sHiveal & hiveal) {

	QPSvcDnaAlignmentRemapper alremap(*this, formValue("parent_proc_ids"),
			formValue("mutualAligmentID"), REALFNLNM, POPULSVCNAME,
			hiveal.dimAl(), 100000);
	alremap.setObjDestination(objs[0].Id());
	return alremap.launch(*user, grpId);
}

idx DnaPopul::OnExecute(idx req) {
	// _/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
	// _/
	// _/ prepare the Alignment file
	// _/
	// _/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
	sHiveId alignerID(formValue("parent_proc_ids"));
	sUsrObj aligner(*user, alignerID);
	sStr path;
	if (aligner.Id()) {
		aligner.getFilePathname00(path,
				"alignment.hiveal" _ "alignment.vioal" __);
	} else {
		reqSetInfo(req, eQPInfoLevel_Error,
				"Object not found or access denied %s\n", alignerID.print());
		reqSetStatus(req, eQPReqStatus_ProgError);
		return 0; // error
	}
	std::auto_ptr<sHiveal> hiveal(new sHiveal(user, path));
	if (!hiveal.get()) {
		reqSetInfo(req, eQPInfoLevel_Error,
				"Object could not be constructed %s\n", alignerID.print());
		reqSetStatus(req, eQPReqStatus_ProgError);
		return 0;
	}
	if (!hiveal->isok() || hiveal->dimAl() == 0) {
		reqSetInfo(req, eQPInfoLevel_Error,
				"Alignment file '%s' is missing or corrupted\n",
				path ? path.ptr() : "unspecified");
		reqSetStatus(req, eQPReqStatus_ProgError);
		return 0;
	}

	// _/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
	// _/
	// _/ load the subject and query sequences
	// _/
	// _/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
	sStr subject;
	aligner.propGet00("subject", &subject, ";");
	sHiveseq Sub(sQPride::user, subject.ptr(), hiveal->getSubMode());
	if (Sub.dim() == 0) {
		reqSetInfo(req, eQPInfoLevel_Error,
				"Reference '%s' sequences are missing or corrupted\n",
				subject ? subject.ptr() : "unspecified");
		reqSetStatus(req, eQPReqStatus_ProgError);
		return 0; // error
	}
	hiveal->Sub = &Sub;

	sStr query;
	aligner.propGet00("query", &query, ";");

	sHiveseq Qry(sQPride::user, query, hiveal->getQryMode());
	if (Qry.dimRefs() > 5)
		Qry.reindex();
	if (Qry.dim() == 0) {
		reqSetInfo(req, eQPInfoLevel_Error,
				"Query '%s' sequences are missing or corrupted\n",
				query ? query.ptr() : "unspecified");
		reqSetStatus(req, eQPReqStatus_ProgError);
		return 0;
	}
	hiveal->Qry = &Qry;

	PERF_START("SORTING ALIGNMENTS");
	// _/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
	// _/
	// _/ fetch all alignments and sorted index of the alignments
	// _/
	// _/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

	if (Sub.dim() > 1) {
		sHiveId mutualAligmentID(formValue("mutualAligmentID"));
		sUsrObj mutAligner(*user, mutualAligmentID);
		if (!mutAligner.Id()) {
			reqSetInfo(reqId, eQPInfoLevel_Error,
					"Multiple alignment object is required.");
			reqSetStatus(reqId, eQPReqStatus_ProgError);
			return 0;
		}
		const char * mutQry = mutAligner.propGet00("query");
		const char * alSub = aligner.propGet00("subject");
		if (sString::size00(mutQry) != sString::size00(alSub)
				|| sString::compareNUntil(mutQry, alSub,
						sString::size00(mutQry), 0, false) == 0) {
			reqSetInfo(reqId, eQPInfoLevel_Error,
					"You must provide a multiple alignment of the subjects used in the alignment step");
			reqSetStatus(reqId, eQPReqStatus_SysError);
			return 0;
		}
		std::auto_ptr<sHiveal> hiveal2(new sHiveal(user));
		sStr reMappedAlignments;
		objs[0].getFilePathname(reMappedAlignments, REALFNLNM);
		hiveal2->parse(reMappedAlignments.ptr());

		if (reMappedAlignments.ptr() && sFile::exists(reMappedAlignments)
				&& hiveal2->isok()) {
			hiveal.reset(hiveal2.release());
			hiveal->Qry = &Qry;
			hiveal->Sub = &Sub;
		} else {
			if (launchRemappers(*hiveal)) {
				reqSetProgress(req, 0, 100);
				reqSetStatus(reqId, eQPReqStatus_Done);
			} else {
				reqSetInfo(reqId, eQPInfoLevel_Error,
						"Failed to launch alignment remapper");
				reqSetStatus(reqId, eQPReqStatus_ProgError);
			}
			return 0;
		}

	}

	idx alCnt = hiveal->dimAl();
	sBioseqAlignment::Al * hdrs = hiveal->getAl(0);

	sVec<idx> alL(sMex::fExactSize);
	alL.add(alCnt);
	idx * alList = alL.ptr(0);

	sVec<idx> sortingOnCoordindateArray(sMex::fExactSize);
	sortingOnCoordindateArray.resize(alCnt);

	idx rd = 1000000;
	sVec<idx> sumVec;
	sumVec.add((1 + ((alCnt - 1) / rd)));
	sumVec.set(0);

	for (idx iAl = 0; iAl < alCnt; ++iAl) {
		idx isV = (iAl / rd);
		hdrs = hiveal->getAl(iAl);
		idx * mS = hiveal->getMatch(iAl);
		sumVec[isV] += hdrs->lenAlign();
		sortingOnCoordindateArray[iAl] = hdrs->subStart() + mS[0];
	}
	real r_avLen = 0;
	for (idx iV = 0; iV < sumVec.dim(); ++iV) {
		r_avLen += sumVec[iV] / alCnt;
	}
	idx avLen = r_avLen;

	sSort::sort(alCnt, sortingOnCoordindateArray.ptr(), alList);
	sortingOnCoordindateArray.destroy();

	PERF_END();

	// _/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
	// _/
	// _/ Initialize frames,consensus vectors, fgprs and alignment clone Ids
	// _/ and windowSize, first position, cutoff for clone calling
	// _/
	// _/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

	sStr bufPath, warningPath00;
	const char * dstPath = bufPath.ptr();

	sVec<int32_t> alClonIds(sMex::fExactSize);
	alClonIds.add(hiveal->Qry->dim());
	alClonIds.set(-1);
	idx /*ascIDs = 0,*/POS = 0, relativeTo = 0;    // stepNum = 0 ;

	bool toEnd = false;
	progress100End = 100;

	real cutoff = formRValue("bifThreshold", 5) / 100;

	reqSetStatus(req, eQPReqStatus_Running);

	bool is_generate_tree = formBoolValue("generateTrees", false);
	// _/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
	// _/
	// _/ loop through all the alignments to produce trees
	// _/
	// _/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
	bool is_stepped_trees = !formBoolValue("InPopulationFrameTrees", false);
	if (is_stepped_trees && is_generate_tree) {
		PERF_START("TREE PLANTING");
		udx tree_frame_size = formUValue("TreeFrameSize", 30);
		udx tree_step_size = formUValue("TreeStep", 30);
		real tree_read_overlap = formRValue("TreeReadOverlap", 1);
		udx tree_size = formUValue("TreeSize", 100);

		udx pool_size = 10 * tree_size;
		sAlTree treeFixedNodes(pool_size, hiveal.get(), alList);
		treeFixedNodes.setFixedMode(0, tree_frame_size, tree_step_size,
				tree_read_overlap, avLen);
		sVec<sAlTree::alNodes> pickedAl;
		pickedAl.empty();
		while (treeFixedNodes.fillRange(pickedAl)) {
			bufPath.add0();
			dstPath = reqAddFile(bufPath, "stepped_result_%" DEC ".tre",
					treeFixedNodes.cur_treeId);
			if (dstPath) {
				sFil newickFile(dstPath);
				if (!pickedAl.dim()
						|| !sAlTree::createNJTree(pickedAl, hiveal.get(),
								newickFile, treeFixedNodes.getStart(),
								treeFixedNodes.getEnd(), tree_size, true)) {
					sFile::remove(dstPath);
				}
			} else {
				warningPath00.printf("stepped_result_%" DEC ".tre",
						treeFixedNodes.cur_treeId);
				warningPath00.add0();
			}
			pickedAl.empty();
		}

		bufPath.add0();
		dstPath = reqAddFile(bufPath, "treesList.jsn");
		if (dstPath) {
			sFil trListFile(dstPath);
			treeFixedNodes.printTreeListJSON(trListFile);
		} else {
			warningPath00.printf("treesList.jsn");
			warningPath00.add0();
		} PERF_END();
	}

	// _/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
	// _/
	// _/ loop through all the alignments
	// _/
	// _/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
	hdrs = hiveal->getAl(0);
	idx start_width = hdrs->lenAlign();

	bool reference_simil = formBoolValue("reference_simil", true);

	bool is_sort_clones = formBoolValue("sortClones", false);

	sVec<idx> statHits, sortedSubjInd;
	if (reference_simil) {
		// _/ Get sorted indexes of subjects based on the alCount. So to output
		// _/ similarities in the same order the appear from alCount.
		sVec<sBioal::Stat> stat;
		hiveal->countAlignmentSummaryBySubject(stat);
		statHits.resize(hiveal->Sub->dim());
		sortedSubjInd.add(hiveal->Sub->dim());
		for (idx is = 0; is < hiveal->dimSub(); ++is) {
			statHits[is] = stat[is + 1].foundRpt;
		}
		sSort::sort(statHits.dim(), statHits.ptr(), sortedSubjInd.ptr());
	}

	sFrameConsensus cloneCons(hiveal->Sub->dim());

	real gap_thr = formRValue("gap_thr", 50) / 100;
	sPopul frame(hiveal.get(), &cloneCons, hiveal->Sub->dim(), 5, alList,
			cutoff, gap_thr);
	frame.fuzzy_threshold = formRValue("fuzzy_threshold", 1) / 100;
	frame._quality_threshold = formIValue("qua_threshold", 0);
	frame.minOverlap = formIValue("minOverlap", 100);
	frame.minBifCov = formIValue("bifMinCov", 0);
	frame.cloneContinuity = formBoolValue("cloneContinuity", false);
	frame.fuzzy_threshold *= cutoff;
	frame.pairEnd_mode =
			formBoolValue("pairedEnd", false) ?
					constructPairMap(frame._pairVec, frame._pairInd,
							hiveal.get(), query.ptr()) :
					0;
	frame.getBifurcationPvalue = formBoolValue("getBifurcationPvalue", false);
	frame.getBifurcationPvalueThrshld = formRValue(
			"getBifurcationPvalueThrshld", 1);
	frame.getBifurcationBayes = formBoolValue("getBifurcationBayes", false);
	frame.getBifurcationBayesThrshld = formRValue("getBifurcationBayesThrshld",
			90) / 100;
	frame.bifurcationDecisionFlag = formIValue("bifurcationTrigger");
	frame.noise_tolerance = formRValue("noise_tolerance", 3) / 100;
	frame.progress_CallbackFunction = sQPrideProc::reqProgressStatic;
	frame.progress_CallbackParam = (void *) this;

	frame.alcl_dist_m = static_cast<sPopul::EAlCloneDistance>(formIValue(
			"distance_measure", static_cast<idx>(sPopul::eEuclidean)));

	frame.InitCoordinates(0, start_width);

	idx frmEND = hiveal->getAl(frame.getAlindex(alCnt - 1))->subStart();
	bool addedAls = false;

	PERF_START("CLONING");
	while (frame.inRange() && !toEnd) {
		// _/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
		// _/
		// _/ Grab the next set of alignments and check if it is the last one
		// _/
		// _/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
		toEnd = frame.getCloneBoundaries(POS, frame.win_size()); //alignments that their subjects start is in range [Start,end) + the fist that start from end;

		// _/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
		// _/
		// _/ Fill windows of clones and bifurcate. Repeat until there is no bifurcation
		// _/
		// _/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
		frame.ensureAliveClone(POS);

		logOut(eQPLogType_Debug,
				"current position %" DEC " ---- Alignment index range [%" DEC ",%" DEC ")\n",
				POS, frame._iAlStart, frame._iAlEnd);

		addedAls = frame.fill(0, &relativeTo);

		if (addedAls) {
			while (frame.bifurcate()) {
				if (!reqProgress(frame.cntAlive(), POS, frmEND))
					return 0;
				frame.fill(1, &relativeTo);
				frame.getBifurcationStats();
			}

//        // _/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
//        // _/
//        // _/ set basecalls for all clones
//        // _/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
			frame.setCalls();

			// _/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
			// _/
			// _/ Merge clones that are similar enough or kill clones
			// _/ with low or zero coverage
			// _/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
			frame.merge();
			frame.kill();
		}

		// _/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
		// _/
		// _/ Decide how big will the next step be.
		// _/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
		idx step = frame.computeStep();

		// _/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
		// _/
		// _/ Extract the consensus of all the clones and define how many positions
		// _/ will be the next step.
		// _/
		// _/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
		frame.extractCons(step, toEnd);

		if (!toEnd) {
			// _/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
			// _/
			// _/ Refill frames up to the position of the next step.
			// _/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
			frame.shift(POS, step);

			// _/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
			// _/
			// _/ Set the next position of the subjects and update progress
			// _/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
			POS += step;
		}
		if (!reqProgress(frame.cntAlive(), POS, frmEND))
			return 0;
	} PERF_END();
	alL.cut(0);
	// _/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
	// _/
	// _/ Output the results
	// _/
	// _/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
	PERF_START("ANALYSIS");

	cloneCons.progress_CallbackFunction = sQPrideProc::reqProgressStatic;
	cloneCons.progress_CallbackParam = (void *) this;

	sVec<idx> cloneIndex;
	sStr erlog;
	bool doBreak = false;
	if (!cloneCons.cleanOutClones(cloneIndex)) {
		erlog.printf("Failed during post processing validation");
		doBreak = true;
	}
	if (!cloneCons.mergeAll()) {
		erlog.printf("Failed to validate merging and bifurcating events");
		doBreak = true;
	}
	if (!cloneCons.trimClones()) {
		erlog.printf("Failed to validate contig boundaries");
		doBreak = true;
	}
	if (!doBreak) {
		cloneCons.setCloneSupport();
		sBioseqpopul::cloneStats clSum;

		cloneCons.getCloneStats(&clSum);

		if (is_sort_clones) {
			cloneCons.reorderClones();
		}

		sVec<idx> ArrSortClon;
		ArrSortClon.add(cloneCons.dim());
		ArrSortClon.set(0);
		idx * arrSortClon = ArrSortClon.ptr(0);
		sBioseqpopul::ParamsCloneSorter paramSort;
		paramSort.flags |= sBioseqpopul::clSortByMaxCov;
		sSort::sortSimpleCallback(sBioseqpopul::sortClonesComparator,
				&paramSort, cloneCons.dim(), cloneCons.ptr(), arrSortClon);

		bufPath.add0();
		dstPath = reqAddFile(bufPath, "clones.simil");
		if (dstPath) {
			cloneCons.sOutSimilarities(dstPath);
		} else {
			warningPath00.printf("clones.simil");
			warningPath00.add0();
		}

		sStr skp2gps, gps2skp, skp2gpsSpprt, gps2skpSpprt;
		idx failCnt = 0, totalcntfls = 4;
		if (!sQPrideProc::reqAddFile(skp2gps, "skip2gap.map")) {
			warningPath00.printf("skip2gap.map");
			warningPath00.add0();
			++failCnt;
		}
		if (!sQPrideProc::reqAddFile(gps2skp, "gap2skip.map")) {
			warningPath00.printf("gap2skip.map");
			warningPath00.add0();
			++failCnt;
		}
		if (!sQPrideProc::reqAddFile(skp2gpsSpprt, "skip2gapSupprted.map")) {
			warningPath00.printf("skip2gapSupprted.map");
			warningPath00.add0();
			++failCnt;
		}
		if (!sQPrideProc::reqAddFile(gps2skpSpprt, "gap2skipSupprted.map")) {
			warningPath00.printf("gap2skipSupprted.map");
			warningPath00.add0();
			++failCnt;
		}
		if (failCnt < totalcntfls) {
			cloneCons.createFrameMaps(skp2gps.ptr(), gps2skp.ptr(),
					skp2gpsSpprt.ptr(), gps2skpSpprt.ptr());
		}

		cloneCons.setConsensusStat(clSum);

		sViopop viopop;
		bufPath.add0();
		dstPath = reqAddFile(bufPath, "clones.viopop");
		if (dstPath) {
			viopop.Digest(dstPath, &cloneCons, arrSortClon, &clSum,
					reference_simil ? sortedSubjInd.ptr() : 0);
		} else {
			erlog.printf("Failed to add result files");
			doBreak = true;
		}
	}

	PERF_END(); PERF_PRINT();

	if (warningPath00.ptr()) {
		bufPath.cut(0);
		warningPath00.add0();
		for (const char * fn = warningPath00.ptr(); fn;
				fn = sString::next00(fn)) {
			bufPath.printf(",%s", fn);
		}
		logOut(eQPLogType_Warning, "Failed to add %s", bufPath.ptr(1));
		reqSetInfo(req, eQPInfoLevel_Warning, "Failed to add %s",
				bufPath.ptr(1));
	}

	if (doBreak) {
		logOut(eQPLogType_Error, erlog.ptr());
		reqSetInfo(req, eQPInfoLevel_Error, erlog.ptr());
		reqSetStatus(req, eQPReqStatus_ProgError);
	} else {
		reqSetProgress(req, frame.cntAlive(), 100);
		reqSetStatus(req, eQPReqStatus_Done);    // change the status
	}

	return 0;
}

int main(int argc, const char * argv[]) {
	sBioseq::initModule(sBioseq::eACGT);

	sStr tmp;
	sApp::args(argc, argv);

	DnaPopul backend("config=qapp.cfg" __,
			sQPrideProc::QPrideSrvName(&tmp, POPULSVCNAME, argv[0]));
	return (int) backend.run(argc, argv);
}
