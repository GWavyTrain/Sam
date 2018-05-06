#!/bin/bash

touch mergers.dat

echo "# tLb [Gyr], z" >> mergers.dat

cat mergers_026-039.dat >> mergers.dat
cat mergers_040-053.dat >> mergers.dat
cat mergers_054-067.dat >> mergers.dat
cat mergers_068-081.dat >> mergers.dat
cat mergers_082-095.dat >> mergers.dat
cat mergers_096-109.dat >> mergers.dat
cat mergers_110-123.dat >> mergers.dat
cat mergers_124-135.dat >> mergers.dat
