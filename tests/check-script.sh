#!/bin/bash


LIST=test.lst
CHK_DIR="./tests/p-adic-checks"

cleanup() {
	rm -rf $LIST 2> /dev/null
	rm -rf $CHK_DIR/*.diff 2> /dev/null
	rm -rf $CHK_DIR/*-new.txt 2> /dev/null
}

cleanup

ls ./bin > $LIST

for i in `cat $LIST`; do
	./bin/$i > $CHK_DIR/$i-new.txt
done

for i in `cat $LIST`; do
	diff $CHK_DIR/$i-new.txt $CHK_DIR/$i.txt > $CHK_DIR/$i.diff
	if [ $? -ne 0 ]; then
		echo "Difference was detected. $i. See:"
		cat $CHK_DIR/$i.diff
		cleanup
		exit 1
	fi
done

echo "OK"

cleanup

exit 0

