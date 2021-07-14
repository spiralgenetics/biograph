#!/bin/bash

SQWISH=`which sqwish`
if [ -z "$SQWISH" ]; then
	echo "Couldn't find sqwish in your PATH. Try 'sudo npm install -g sqwish'"
	exit 1
fi

UGLIFYJS=`which uglifyjs`
if [ -z "$UGLIFYJS" ]; then
	echo "Couldn't find uglifyjs in your PATH. Try 'sudo npm install -g uglify-js'"
	exit 1
fi

set -e

if [ ! -r 'js/drawchart.js' ]; then
  echo "Please cd to the kmer_quality_report source tree first."
  exit 1
fi

if [ ! -z "$1" ]; then
	echo "Enabling unit tests."
	cp tests/10-tests.html compiled/footer/10-tests.html
	cp tests/qunit.css css/qunit.css
else
	echo "Disabling unit tests. Enable with '$0 test'."
	rm -f compiled/footer/10-tests.html
	rm -f css/qunit.css
fi

JQUERY='jquery-2.1.0.min.js'
#JQUERY='jquery-1.11.0.min.js'

JQUERYUI='jquery-ui-1.10.4.custom.min.js'
BOOTSTRAP='bootstrap.min.js'
UNDERSCORE='underscore.min.js'
D3='d3.v3.min.js'
NVD3='nv.d3.js'
CROSSFILTER='crossfilter.min.js'

echo "Compressing css..."
cat css/*.css > all-css.css
$SQWISH all-css.css
mv all-css.min.css compiled/header/02-css.css
rm -f all-css.css

cd js
echo > all-libraries.js
for x in $JQUERY $JQUERYUI $UNDERSCORE $BOOTSTRAP $D3 $NVD3 $CROSSFILTER; do
	cat $x >> all-libraries.js
	echo >> all-libraries.js
done
cd ..

if [ -z "$QUICK" ]; then
	WARNINGS='warnings.txt'
	echo "Compressing libraries..."
	$UGLIFYJS --screw-ie8 -c -m -o compiled/header/04-header-scripts.js js/all-libraries.js 2> $WARNINGS
	WCOUNT=`grep -c . $WARNINGS`
	if [ $WCOUNT -gt 0 ]; then
		echo "$WCOUNT compilation warnings. See $WARNINGS for details."
	fi
else
	echo "Skipping library compression."
fi

echo "Compressing app code..."
$UGLIFYJS --screw-ie8 -c -m -o compiled/footer/03-app.js js/drawchart.js

echo "Generating kmer_quality_report.html..."
cat compiled/header/* compiled/data.js compiled/footer/* > kmer_quality_report.html

if [ -z "$1" ]; then
	echo "Making deployment files..."
	cat compiled/header/* compiled/footer/* > kmer_quality_report.h
fi
