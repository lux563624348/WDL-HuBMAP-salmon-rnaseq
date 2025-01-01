#dockstore workflow launch --local-entry ../steps/fastqc.wdl --json test.wdl.json
miniwdl run ./pipeline.wdl -i ./test/test.wdl.json
