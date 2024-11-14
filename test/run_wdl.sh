#dockstore workflow launch --local-entry ../steps/fastqc.wdl --json test.wdl.json
miniwdl run ../steps/test.wdl -i test.wdl.json
