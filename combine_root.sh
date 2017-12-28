#/bin/bash
MRB_T="/uboone/app/users/ilker228/v06_26_00/srcs/uboonecode/uboone/AnalysisTree"
find $PWD -name '*.root' >$MRB_T/Final.txt
hadd -f $MRB_TOP/Grid/newilkerdata.root @$MRB_T/Final.txt
