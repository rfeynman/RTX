
allprog=mainq mainp mainga preParmela  lintrack opttrack oqttrack maint longPlot astraDist survey fixOrbit rtxupdate  readTraj qplot
allexe=mainq.exe mainp.exe mainga.exe  preParmela.exe  lintrack.exe opttrack.exe oqttrack.exe maint.exe longPlot.exe astraDist.exe survey.exe qplot.exe
allexe=$(allprog)
all:	$(allprog)



readTraj:	readTraj.cxx
		g++ -o readTraj readTraj.cxx

readGDF:	readGDF.cxx
		g++ -o readGDF readGDF.cxx


rtxupdate:	rtxupdate.cxx 
	g++ rtxupdate.cxx -g -o rtxupdate 

fixOrbit:	fixOrbit.cxx 
	g++ fixOrbit.cxx -g -o fixOrbit 

survey:	survey.cxx 
	g++ survey.cxx -g -o survey 

oqttrack:	oqttrack.cxx pmla.o  dminv.o
	g++ oqttrack.cxx -g -o oqttrack pmla.o dminv.o -I$(HOME)/condorCommonN  $(HOME)/condorCommonN/libcondor.a -lpthread


opttrack:	opttrack.cxx pmla.o  dminv.o
	g++ opttrack.cxx -g -o opttrack pmla.o dminv.o -I$(HOME)/condorCommonN  $(HOME)/condorCommonN/libcondor.a -lpthread

lintrack:	lintrack.cxx pmla.o  dminv.o
	g++ lintrack.cxx -g -o lintrack pmla.o dminv.o 

mainq:	mainq.cxx pmla.o dminv.o
	g++ -o mainq mainq.cxx pmla.o dminv.o -I$(HOME)/condorCommonN $(HOME)/condorCommonN/libcondor.a -lpthread -g

mainp:	mainp.cxx pmla.o dminv.o
	g++ -o mainp mainp.cxx pmla.o dminv.o -I$(HOME)/condorCommonN $(HOME)/condorCommonN/libcondor.a -lpthread -g

mainga:	mainga.cxx pmla.o dminv.o
	g++ -o mainga mainga.cxx pmla.o dminv.o -lpthread -g

gptopt:	gptopt.cxx 
	g++ -o gptopt gptopt.cxx -I$(HOME)/condorCommonN $(HOME)/condorCommonN/libcondor.a -lpthread -g

maint:	maint.cxx pmla.o dminv.o
	g++ -o maint maint.cxx pmla.o dminv.o -I$(HOME)/condorCommonN $(HOME)/condorCommonN/libcondor.a -lpthread -g

mainm:	mainm.cxx  dminv.o
	g++ -o mainm mainm.cxx  dminv.o -I$(HOME)/condorCommon $(HOME)/condorCommon/libcondor.a -g


longPlot:	longPlot.cxx pmla.o dminv.o xmgr_class.o 
	g++ -o longPlot longPlot.cxx pmla.o dminv.o xmgr_class.o -g

astraDist:	astraDist.cxx pmla.o dminv.o xmgr_class.o 
	g++ -o astraDist astraDist.cxx pmla.o dminv.o xmgr_class.o -g

qplot:	qplot.cxx pmla.o dminv.o xmgr_class.o 
	g++ -o qplot qplot.cxx pmla.o dminv.o xmgr_class.o -g


preParmela:	preParmela.cxx sobseq.cxx 
	g++ preParmela.cxx sobseq.cxx -g -o preParmela  -I$(HOME)/condorCommon  $(HOME)/condorCommon/libcondor.a





dminv.o:	dminv.cxx
	g++ dminv.cxx -c -g

poly.o:	poly.cxx poly.hxx
	g++ poly.cxx -c -g


pmla.o:	pmla.cxx pmla.hxx
	g++ pmla.cxx -c -g


xmgr_class.o:	xmgr_class.cxx xmgr_class.hxx
	g++ xmgr_class.cxx -c -g




clean:
	rm -f *.o $(allexe)
