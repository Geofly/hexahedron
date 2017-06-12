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
#pragma once
#ifndef QPSVC_DNAALIGNMENTREMAPPER_hpp
#define QPSVC_DNAALIGNMENTREMAPPER_hpp

#include <qpsvc/qpsvc.hpp>

class QPSvcDnaAlignmentRemapper: public sQPSvc
{
        typedef sQPSvc TParent;

    public:

        QPSvcDnaAlignmentRemapper(sQPride & qp, const char *alId, const sHiveId &msaId, const char * resultsFilename, const char * serviceToBeLaunched, udx alCnt , udx chnkSize)
            : TParent(qp), m_totSize(alCnt), m_maxChunkSize(chnkSize)
        {
            setAlignmentId(alId);
            setMultipleAlignmentId(msaId);
            setResultsFilename(resultsFilename);
            setServiceToBeLaunched(serviceToBeLaunched);
        }

        virtual ~QPSvcDnaAlignmentRemapper()
        {
        }

        void setAlignmentId(const sHiveId & objId)
        {
            setVar("alignment", "%s", objId.print());
        }
        void setMultipleAlignmentId(const sHiveId & qryId)
        {
            setVar("multipleAlignment", "%s", qryId.print());
        }
        virtual udx split() const
        {
            if(!m_maxChunkSize)
                return 1;
            return 1+((m_totSize-1)/m_maxChunkSize);
        }
        void setResultsFilename(const char * file)
        {
            setVar("resultsFileName", "%s", file);
        }
        void setServiceToBeLaunched(const char * serviceToBeLaunched)
        {
            setVar("serviceToBeLaunched", "%s", serviceToBeLaunched);
        }
        void setObjDestination(const sHiveId & objId)
        {
            setVar("obj", "%s", objId.print());
        }
        virtual const char* getSvcName() const
        {
            return "dna-alignment-remapper";
        }

        udx m_totSize;
        udx m_maxChunkSize;
};

#endif
