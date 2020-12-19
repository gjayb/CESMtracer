#!/bin/csh -f

#------------------------------------------------------------------------------------
# For now, set streams manually. You must only set as many streams as are declared
#  in the tavg_nml section. For example, if there are three streams:
#  @ s1 = $my_stream
#  @ s2 = $s1 + 1
#  @ s3 = $s2 + 1
#------------------------------------------------------------------------------------

@ my_stream = $1
if ($my_stream < 1) then
   echo invalid my_stream number $my_stream
   exit 5
endif

@ s1 = 1   # use base-model stream 1
@ s2 = 2   # use base-model stream 2

cat >! $CASEROOT/Buildconf/popconf/bipit_tavg_contents << EOF
$s1  TAUML
$s1  NUTRI
$s1  PARTI
$s1  J_NUTRI 
EOF

#@ s2 = $s1 + 1   # use base-model stream 2

#cat >! $CASEROOT/Buildconf/popconf/bipit_tavg_contents << EOF
#$s1  TAUML
#$s1  NUTRI
#EOF
#-------------------------------------------------------------------------------------
# Add optional tracer budget terms
#-------------------------------------------------------------------------------------
if ($OCN_TAVG_TRACER_BUDGET == TRUE) then
cat >> $CASEROOT/Buildconf/popconf/bipit_tavg_contents << EOF
$s1  BIPIT_RESET_TEND
$s1  J_NUTRI
$s1  Jint_100m_NUTRI
$s1  tend_zint_100m_NUTRI
$s1  UE_NUTRI
$s1  VN_NUTRI
$s1  WT_NUTRI
$s1  HDIFE_NUTRI
$s1  HDIFN_NUTRI
$s1  HDIFB_NUTRI
$s1  DIA_IMPVF_NUTRI
$s1  KPP_SRC_NUTRI
$s1  J_PARTI
$s1  Jint_100m_PARTI
$s1  tend_zint_100m_PARTI
$s1  UE_PARTI
$s1  VN_PARTI
$s1  WT_PARTI
$s1  HDIFE_PARTI
$s1  HDIFN_PARTI
$s1  HDIFB_PARTI
$s1  DIA_IMPVF_PARTI
$s1  KPP_SRC_PARTI
$s1  J_TAUML
$s1  Jint_100m_TAUML
$s1  tend_zint_100m_TAUML
$s1  UE_TAUML
$s1  VN_TAUML
$s1  WT_TAUML
$s1  HDIFE_TAUML
$s1  HDIFN_TAUML
$s1  HDIFB_TAUML
$s1  DIA_IMPVF_TAUML
$s1  KPP_SRC_TAUML
EOF
endif

#  disable the following until they are computed correctly
#  BIPIT_SQR 
#  UE_BIPIT
#  VN_BIPIT
#  WT_BIPIT
#  ADV_BIPIT
#  J_BIPIT
#  Jint_BIPIT
#  STF_BIPIT
#  RESID_BIPIT
#  FvPER_BIPIT
#  FvICE_BIPIT

