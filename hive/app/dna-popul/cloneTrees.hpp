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
#ifndef cloneTrees_hpp
#define cloneTrees_hpp

#include <slib/utils.hpp>
#include <ssci/math/clust/clust2.hpp>

#include <ssci/bio.hpp>
#include <violin/violin.hpp>

#define DEBUGGING_KILLING
#define DEBUGGING_MERGING

#define RAND1 (real)rand()/RAND_MAX

class sAlTree {

    public:
        struct cloneNode {
            idx posStart;
            idx posEnd;
            idx iAl;
            sVec<idx> treeList;
        };
        struct treeStat {
            idx id, start, end, dim,startAl, endAl;
        };
        struct alNodes {
            idx Al;
            idx rpt;
        };
        sVec<cloneNode> readDic;
        sVec<cloneNode> readsList;
        sVec<treeStat> treesList;
        idx cur_treeId;
        idx _size,_space;
        real _overlap;
        real _forceRefreshPortion;

    private:
        idx _startAl,_endAl;
        idx _start,_end;
        idx _fStart,_fEnd; //theoretical first end (start+readSize) and last start (end- readSize)
        idx _oStart, _oEnd;
        idx _iStart, _iEnd;
        idx * _alList; //sorted array
        idx _step, _width, _curWidth, _readSize; //_width given by user _curWidth considering depth and coverage
        real _readOffset;
        sHiveal * _hiveal;
        bool _isFixedMode;
        bool _isRandom;

        // Setting Alingment Indices
        void _setEndAl(idx endAl){ //exclusive
            _endAl = endAl;
        }
        void _setStartAl(idx startAl){ //inclusive
            _startAl = startAl;
        }

        // Validation checks
        bool _isNodeInOverlap(cloneNode * node) {
            return _isInOverlap(node->posStart, node->posEnd);
        }

        bool _isNodeInFrame (cloneNode * node) {
            return _isInFrame(node->posStart, node->posEnd);
        }

        bool _isNodeValid(cloneNode * node) {
            return _isValid(node->posStart, node->posEnd);
        }

        bool _isInOverlap (idx st, idx end) {

            bool res = ( st <= _oStart ) && ( end >= _oEnd );
            return res;
        }

        bool _isInFrame (idx st, idx end) {
            idx t_start = (_iStart>=0)?_iStart:0;
            idx t_end = _iEnd ;
            bool res = ( end > t_start ) && (st <= t_start);
            if ( !res ) {
                res = ( st < t_end ) && ( end >= t_end );
            }
            return res;
        }

        bool _isValid (idx st, idx end) {
            if( !_overlap ) {
                return _isInFrame(st,end);
            }
            else {
                return _isInFrame(st,end) && _isInOverlap(st,end);
            }
        }


        //Setting basic variables
        void _setWidth (idx width) {
            _width = width;
        }

        void _setCord ( idx start ) {
            _fStart = start;
            _start = start;
            _fEnd = start + _width;
            _end = _fEnd;

        }

        void _setStep (idx step ) {
            _step = step;
        }

        void _setReadSize (idx readSize) {
            _readSize = readSize;
        }


        // Getting overlap info
        idx _getOverlapStart () {
            idx t_width = _getOverlapSize ();
            return _start + (t_width - t_width*_overlap)/2;
        }
        idx _getOverlapSize () {
            idx width = _width;
            if( _width > _curWidth && _curWidth )
                width = _curWidth;
            return width * _overlap;
        }
        idx _getOverlapEnd () {
            return _getOverlapStart() + _getOverlapSize();
        }
        idx _getAl (idx iAl) {
            return  _alList?_alList[iAl]:iAl;
        }


        //Setting overlap variables
        void _setOverlap (real overlap) {
            if( overlap > 1 ) {
                overlap = 1;
            }
            else if ( overlap < 0 ) {
                overlap = 0;
            }
            _overlap = overlap;
        }
        void _setOverlapCord () {
            if( !_overlap ) {
                _oStart = _start;
                _oEnd = _end;
            }
            else {
                _oStart = _getOverlapStart();
                _oEnd = _getOverlapEnd();
            }
        }

        void _setInputOffset () {
            if( !_overlap ) {
                _iStart = _start;
                _iEnd = _end;
            }
            else {
                idx dif = _readSize - _getOverlapSize();
                if( dif < 0 ) dif = 0;
                _iStart = _oStart - dif;
                _iEnd = _start+1;
            }
        }


        idx _generateRandVector(idx iStart, idx iEnd,idx cnt, sVec<idx> & out);
        idx _getPosEnd(idx iAl);
        idx _getPosStart(idx iAl);
        idx  _getNextAlInds( bool reStart = false );
        idx _refreshPoolToPos();
        void _updateStartsOnAl();
        void _updateEndsOnAl();
        void _moveFrameInds();
        void _addTree() ;
        void _dic2VecIDs(sVec<alNodes> & qryIds) ;


    public:

        sAlTree ( idx size,sHiveal * hiveal, idx * alList=0){
            _size = size;
            _space = _size;
            _alList = alList;
            _hiveal = hiveal;
            init();
        }
        void init() {
            _startAl = 0;_endAl = 0;
            _start = -1;_end = -1;
            _oStart = -1;_oEnd = -1;
            _iStart = -1;_iEnd = -1;
            _fStart = -1; _fEnd = -1;
            _isFixedMode = false;
            _isRandom = false;
            readDic.empty();
            treesList.empty();
            readsList.empty();
            cur_treeId=0;
            _space =  _size ;
            _step = 0;
            _overlap = 0;
            _width = 0;
            _curWidth = 0;
            _forceRefreshPortion = 0.5;
        }
        void setFixedMode (idx start, idx width, idx step, real overlap, idx readSize) {
            _isFixedMode = true;
            _setWidth( width );
            _setStep( step );
            _setOverlap ( overlap );
            _setReadSize ( readSize );

            _setCord ( start-step );
            _setOverlapCord();
            _setInputOffset ();
        }

        idx getStart() {
            return _start;
        }
        idx getEnd() {
            return _end;
        }

        idx appendToReadsList(cloneNode * node);

        idx fillRange(sVec<alNodes> & qryIds, idx * startAl = 0 , idx * endAl = 0 );

        void printTreeListJSON( sStr & buf);

        bool isFull() {
            return !_space;
        }
        bool isRandom() {
            return _isRandom;
        }
        bool isFixedMode() {
            return _isFixedMode;
        }

        static idx createNJTree(sVec<sAlTree::alNodes> & mappedAl,sBioal * bioal, sStr &buf, idx frameStart, idx frameEnd, idx treeSize = 0, bool fixedTreeSize = false);
        static idx calculatePairwiseDistance(sBioal * bioal, idx * ai1, idx * ai2, idx frameStart, idx frameEnd, bool countNonOverlappingRegions = true);
        static idx nodePrintfCallback(sStr &out, sHierarchicalClustering &clust, idx x, void *param)
        {
            sStr * idList = static_cast<sStr*>(param);
            if (x >= 0 && x < sString::cnt00( idList->ptr() ))
                out.printf("%s", sString::next00(idList->ptr(),x) );
            return 0;
        }
};
#endif //cloneTrees_hpp
