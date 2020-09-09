#!/bin/bash
scp TimtheTyrant@login05.osgconnect.net:pde-sim-transfer/Log.tar.gz ../data
cd ../data/
tar -xzf Log.tar.gz
rm Log.tar.gz