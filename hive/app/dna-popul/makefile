#/*
# * $Id$ ::2B6A7A!
# * Copyright (c) 2005 Dr. Vahan Simonyan and Dr. Raja Mazumder.
# * This software is protected by U.S. Copyright Law and International
# * Treaties. Unauthorized use, duplication, reverse engineering, any
# * form of redistribution, use in part or as a whole, other than by
# * prior, express, written and signed agreement is subject to penalties.
# * If you have received this file in error, please notify copyright
# * holder and destroy this and any other copies. All rights reserved.
# */

PROGRAMS += dna-popul

OBJECT_TYPES = svc-popul

PROGRAMS_LIBS = -lviolib -lqpsvc $(VLIB_LIB)
PROGRAMS_XLIBS = $(MYSQL_LIB)

dna-popul_SRCS = dna-popul.cpp bioclones.cpp cloneTrees.cpp

ifeq ($(devdir),)
    devdir = $(abspath ../../..)
endif
include $(devdir)/build/makefile.meta

PNOHY: dist
dist::
	$(call dist-service-register,dna-popul,HIVE-HexaHedron QSR Analysis tool,/Sequencing/)
	$(call dist-app,dna-popul,dna-popul)
