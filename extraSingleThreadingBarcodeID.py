def singleThreadedReading():
	for j in windowList[0:4]:
		print(j)
		tmpDict = {}
		aboveMinCutoffDict = {}
		indFastq = "readsOn."+j+".fastq"
		ecFastq = errorCorrectDir+"/"+j+"/"+"readsOn."+j+".ec.fastq"
		if os.path.isfile(ecFastq):
			pass
		else:
			print "ln -s " + indFastq + " " + ecFastq
			subprocess.call("ln -s " + indFastq + " " + ecFastq, shell=True)
	
	
		fastqDict = SeqIO.to_dict(SeqIO.parse(ecFastq,"fastq"))

		counter = 0
		for key in fastqDict:
			counter += 1 
			seqK = fastqDict[key].seq
			if seqK in tmpDict:
				tmpDict[seqK] = tmpDict[seqK] + 1
			else:
				tmpDict[seqK] = 1
	
		sumReads = len(fastqDict)
		sumReadFilteredReads = sumReads
		# tmpDict is BC sequence with times seen per window
		if bool(tmpDict):
			print(tmpDict)
			coverageList = list(tmpDict.values())

			percLimit = findSumThreshold(coverageList, percentBCReads)
			print(percLimit)

			tmp2Dict = {}
			for k, v in tmpDict.iteritems():
				if v >= percLimit:
					tmp2Dict[k] = v

			for tmp2Key in tmp2Dict:
				if tmp2Dict[tmp2Key] >= minBarcodeCounts:
					aboveMinCutoffDict[tmp2Key] = tmp2Dict[tmp2Key]
		
			for tmp3Key in aboveMinCutoffDict:
				tmp3List = [j, aboveMinCutoffDict[tmp3Key]]
				if tmp3Key in bcDictSeqsAsKeys:
					bcDictSeqsAsKeys[tmp3Key].append(tmp3List)
				else:
					bcDictSeqsAsKeys[tmp3Key] = list()
					bcDictSeqsAsKeys[tmp3Key].append(tmp3List)
		else:
			continue