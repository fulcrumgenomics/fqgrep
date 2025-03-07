#!/bin/bash

set -euo pipefail

echo -e "<!-- start usage -->\n\`\`\`console\n" > usage.txt;
# the `sed` is needed since the # of cores detected may differ across environments
./target/debug/fqgrep --help | ./target/debug/fqgrep --help | sed -e '\/[default: [0-9]*\]/ s/\[default: [0-9]*\]/[default: 10]/' >> usage.txt;
echo -e "\`\`\`\n<!-- end usage -->" >> usage.txt;
sed -e '/<!-- start usage -->/,/<!-- end usage -->/!b' -e '/<!-- end usage -->/!d;r usage.txt' -e 'd' README.md > README.md.new;
mv README.md.new README.md;
rm usage.txt;
