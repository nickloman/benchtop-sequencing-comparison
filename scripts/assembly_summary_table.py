import sys
from runutils import read_assemblies, read_run_details, hashit

assemblies = read_assemblies(sys.argv[1])
summaries = hashit(read_run_details(sys.argv[2]), 'Name')

#Name   NumContigs      NumRefReplicons NumAssemblyBases        NumReferenceBases       NumLCBs DCJ_Distance    NumDCJBlocks    NumSNPs NumMisCalled    NumUnCalled     NumGapsRef      NumGapsAssembly TotalBasesMissed        PercBasesMissed ExtraBases      PercExtraBases  MissingChromosomes      ExtraContigs    NumSharedBoundaries     NumInterLcbBoundaries   BrokenCDS       IntactCDS       ContigN50       ContigN90       MinContigLength MaxContigLength AA      AC      AG      AT      CA      CC      CG      CT      GA      GC      GG      GT      TA      TC      TG      TT

fields = ['NumContigs', 'NumAssemblyBases', 'MaxContigLength', 'ContigN50', 'NumLCBs', 'NumGapsRef', 'NumGapsAssembly', 'PercBasesMissed']

print "Sample" + "\t" + "Assembler" + "\t" + "\t".join(fields)

for a in assemblies:
	try:
		s = summaries[a['Name'] + '.fas']
	except:
		s = summaries[a['Name']]
	print "%s\t%s\t" % (a['Desc'], a['AssemblySoftware']) ,
	print "\t".join([s[f] for f in fields])

