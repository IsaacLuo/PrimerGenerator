import Bio.SeqIO
import json

def setParameter():
    try:
        para = json.load(open('parameter.json'))

    except:
        print 'parameter file is damage'

    for p in para:
        print "%s: %s"%(p,para[p])

    yn = raw_input("continue?(y/n)")

    if yn.lower() == 'y':
        print 'start'
        return para
    else:
        print 'quit'
        return None

def getTM2(seq):
    return 2*(seq.count('A') + seq.count('T')) + 4*(seq.count('C') + seq.count('G'))

def getTM(seq):
    ## Tm= 64.9 +41*(yG+zC-16.4)/(wA+xT+yG+zC)
    return 64.9 + 41*(seq.count('G')+seq.count('C')-16.4)/len(seq)

def checkEnd(seq):
    return seq[-1] == 'C' or seq[-1] == 'G'

def getGC(seq):
    return 100*(seq.count('C') + seq.count('G')) / (len(seq)*1.0)

def generateIntegerRange(lowV,highV,optV):
    if highV<lowV or optV < lowV or optV > highV:
        return []
    re = [optV]
    maxOffset = max(optV-lowV,highV-optV)
    for offset in range(1,maxOffset+1):
        if optV + offset <= highV:
            re.append(optV+offset)
        if optV - offset >= lowV:
            re.append(optV-offset)
    return re

def designSequencingPrimers(fileName, para):
    f = open(fileName)
    seqRecord = Bio.SeqIO.read(f,'fasta')
    f.close()
    
    seq = seqRecord.seq.upper()

    primerLengthRange = generateIntegerRange(para['primerLengthRange'][0],para['primerLengthRange'][1],para['primerLengthOptimize'])
    startPos = para['startPosition']
    endPos = para['endPosition'] if para['endPosition'] > 0 and para['endPosition'] > para['startPosition'] else len(seq)
    tmRange = para['tmRange']
    gcRange = para['gcRange']

    segmentRange = generateIntegerRange(para['segmentRange'][0]-para['overlapRange'][1],para['segmentRange'][1]-para['overlapRange'][0],para['segmentOptimize']-para['overlapOptimize'])
    
    maxOffsetTry = len(segmentRange)
    
    
    primerPoses = []
    print 'len:%d'%len(seq)
    s0 = startPos - (para['segmentOptimize']-para['overlapOptimize'])
    searchBack = 0;
    while s0 < endPos - (para['segmentOptimize']-para['overlapOptimize']):
        for segmentLength in segmentRange:
            found = False
            s = s0+segmentLength
            print "%d+%d = %d"%(s0,segmentLength,s)
            for l in primerLengthRange:
                primer = seq[s:s + l]
                if checkEnd(primer):
                    tm = getTM(primer)
                    gc = getGC(primer)
                    if tm>=tmRange[0] and tm<=tmRange[-1] and gc>=gcRange[0] and gc<=gcRange[-1]:
                        print 'found'
                        primerPoses.append([s,s+l])
                        pos = [s,s+l]
                        s0 = s
                        searchBack = 0;
                        found = True
                        break
            if found==True:
                break
        else:
            print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            break
                        
    #print primerPoses
    if (len(primerPoses)>0):
        with open(fileName+'.csv','w') as fcsv:
            fcsv.write('start,end,offset,length,tm,gc,sequence\n')
            lastStart = primerPoses[0][0]
            for pos in primerPoses:
                print '[%d-%d][OFF:%d][LEN:%d][TM:%d][GC:%.1f%%] %s'%(pos[0],pos[1],pos[0]-lastStart,pos[1]-pos[0],getTM(seq[pos[0]:pos[1]]),getGC(seq[pos[0]:pos[1]]),seq[pos[0]:pos[1]])
                fcsv.write('%d,%d,%d,%d,%d,%.1f%%,%s\n'%(pos[0],pos[1],pos[0]-lastStart,pos[1]-pos[0],getTM(seq[pos[0]:pos[1]]),getGC(seq[pos[0]:pos[1]]),seq[pos[0]:pos[1]]))
                lastStart = pos[0]
                
            
    f.close()
 


parameter = setParameter()
if parameter != None:
    for name in parameter['fileName']:
        designSequencingPrimers(name,parameter)

           
