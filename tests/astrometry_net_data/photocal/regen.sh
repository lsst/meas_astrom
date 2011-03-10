#!/usr/bin/env bash

ORIG=index-2033.orig.fits
NEW=index-2033.fits
FITSCOPY=$ASTROMETRY_NET_DIR/bin/fitscopy
BUILDIDX=$ASTROMETRY_NET_DIR/bin/build-index

$FITSCOPY ${ORIG}+13"[col mag;id;mag_err(e)=0.1;starnotgal(L)=T;variable(L)=F]" $NEW

