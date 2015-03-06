#!/bin/bash
if [ ! -e .end ];then
    {
	rm evolution.dump
	touch .start
	date +%s > .time
	./tidev-full.out
	date +%s >> .time
	touch .end
	echo "Done."
    }|tee run.log
else
    echo "Already ran."
fi
