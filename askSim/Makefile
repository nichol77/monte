#############################################################################
## Makefile -- New Version of my Makefile that works on both linux
##              and mac os x
## Ryan Nichol <rjn@hep.ucl.ac.uk>
##############################################################################
include Makefile.arch

#Site Specific  Flags
SYSINCLUDES	= -I/sw/include
SYSLIBS         = -L/sw/lib

#Now the bits we're actually compiling
ROOT_LIBRARY = libAskRay.${DllSuf}
LIB_OBJS = InteractionInfo.o NeutrinoEvent.o AskConventions.o AskGeom.o BedmapTable.o WorldModel.o SecondaryPropagator.o interactions.o NeutrinoGenerator.o UsefulPhysics.o UtilityFuncs.o NeutrinoCrossSectionModel.o NeuGenRunSummary.o RFSensor.o RFPulse.o askDict.o 
CLASS_HEADERS = InteractionInfo.h NeutrinoEvent.h AskConventions.h AskGeom.h BedmapTable.h WorldModel.h SecondaryPropagator.h NeutrinoGenerator.h UsefulPhysics.h NeutrinoCrossSectionModel.h UtilityFuncs.h NeuGenRunSummary.h RFSensor.h RFPulse.h

PROGRAMS = generateSecondaries generateNeutrinos

all : $(ROOT_LIBRARY) $(PROGRAMS)

#The programs
generateNeutrinos : $(ROOT_LIBRARY) generateNeutrinos.$(SrcSuf)
	@echo "<**Compiling**> "  
	$(LD)  $(CXXFLAGS) $(LDFLAGS) generateNeutrinos.$(SrcSuf) $(ROOT_LIBRARY) $(LIBS) -o $@

generateSecondaries : $(ROOT_LIBRARY) generateSecondaries.$(SrcSuf)
	@echo "<**Compiling**> "  
	$(LD)  $(CXXFLAGS) $(LDFLAGS) generateSecondaries.$(SrcSuf) $(ROOT_LIBRARY) $(LIBS) -o $@



#The library
$(ROOT_LIBRARY) : $(LIB_OBJS) 
	@echo "Linking $@ ..."
ifeq ($(PLATFORM),macosx)
# We need to make both the .dylib and the .so
		$(LD) $(SOFLAGS)$@ $(LDFLAGS) $^ $(OutPutOpt) $@
ifneq ($(subst $(MACOSX_MINOR),,1234),1234)
ifeq ($(MACOSX_MINOR),4)
		ln -sf $@ $(subst .$(DllSuf),.so,$@)
else
		$(LD) -bundle -undefined $(UNDEFOPT) $(LDFLAGS) $^ \
		   $(OutPutOpt) $(subst .$(DllSuf),.so,$@)
endif
endif
else
	$(LD) $(SOFLAGS) $(LDFLAGS) $(LIBS) $(LIB_OBJS) -o $@
endif


%.$(ObjSuf) : %.$(SrcSuf)
	@echo "<**Compiling**> "$<
	$(CXX) $(CXXFLAGS) -c $< -o  $@

%.$(ObjSuf) : %.C
	@echo "<**Compiling**> "$<
	$(CXX) $(CXXFLAGS) $ -c $< -o  $@


askDict.C: $(CLASS_HEADERS)
	@echo "Generating dictionary ..."
	@ rm -f *Dict* 
	rootcint $@ -c $(CLASS_HEADERS) LinkDef.h

clean:
	@rm -f *Dict*
	@rm -f *.${ObjSuf}
	@rm -f $(LIBRARY)
	@rm -f $(ROOT_LIBRARY)
	@rm -f $(subst .$(DllSuf),.so,$(ROOT_LIBRARY))	
	@rm -f $(TEST)
#############################################################################



