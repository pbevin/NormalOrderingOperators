#!/bin/sh

python3 NOO.py --selftest || exit 1
python3 NOO.py 'ad1s[k1].ad1s[k2].ah1s[p1mu].cd1s[p1m].ad1s[p2m].cd1s[q1].cd1s[q2].ch1s[q3]' > /tmp/result.txt \
  && diff test_ans.txt /tmp/result.txt \
  && echo "OK"
