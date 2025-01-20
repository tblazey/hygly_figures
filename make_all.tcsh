#!/bin/tcsh -f

foreach fig ( `ls -d demo table* figure_* ` )

    cd $fig
    tcsh ./${fig}.tcsh
    if ( $? == 1 ) exit
    cd ../

end