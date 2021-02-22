#!/bin/sh
START=$(date +%s)
# do something
# start your script work here
./Run.bash
ls -R /etc > /tmp/x
rm -f /tmp/x
# your logic ends here
END=$(date +%s)
DIFF=$(( $END - $START ))
echo "It took $DIFF seconds"
