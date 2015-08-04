#!/bin/bash

for i in /var/www/html/phyloprofile/tmpData/*
do
    # file date
    fileDate=`stat "$i" | grep Modify | awk '{print $2, $3}'`
    transDate=`date -d "$fileDate" +%s`

    # 2 hours ago
    ago2h=`date -d '-2 hours' +%s`

    # compare file date
    if [ $transDate -le $ago2h ]; then
	`rm -rf "$i"`
    fi
done


