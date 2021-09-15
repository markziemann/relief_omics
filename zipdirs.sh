#!/bin/bash
for DIR in `ls -d */ | sed 's#/##' ` ; do
  ZIP=$DIR.zip ; zip -r $ZIP $DIR/ &
done
