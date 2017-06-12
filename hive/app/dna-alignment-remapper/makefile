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

PROGRAMS += dna-alignment-remapper

OBJECT_TYPES = svc-alignment-remapper

PROGRAMS_LIBS = -lviolib $(VLIB_LIB)
PROGRAMS_XLIBS = $(MYSQL_LIB)

dna-alignment-remapper_SRCS = dna-alignment-remapper.cpp

ifeq ($(devdir),)
    devdir = $(abspath ../../..)
endif
include $(devdir)/build/makefile.meta

PNOHY: dist
dist::
	$(call dist-service-register,dna-alignment-remapper,HIVE-Alignment Remapper tool,/Sequencing/)
	$(call dist-app,dna-alignment-remapper,dna-alignment-remapper)

