#!/bin/bash

for i in /var/www/html/phyloprofile/tmpData/*
do
    # file date
    fileDate=`stat "$i" | grep Modify | awk '{print $2, $3}'`
    transDate=`date -d "$fileDate" +%s`

    # 6 hours ago
    ago6h=`date -d '-2 hours' +%s`

    # compare file date
    if [ $transDate -le $ago6h ]; then
	`rm -rf "$i"`
    fi
done


