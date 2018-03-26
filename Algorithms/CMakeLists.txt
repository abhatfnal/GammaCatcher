art_make(
          LIB_LIBRARIES
			larcore_Geometry	
		   	uboone_Geometry		
			lardata_Utilities	
			larcore_Geometry_Geometry_service	
			lardataobj_RecoBase
			lardataobj_AnalysisBase	
			canvas_Utilities
                        ${MF_MESSAGELOGGER}
                        ${MF_UTILITIES}
                        ${CETLIB}
                        ${Boost_SYSTEM_LIBRARY}
                        ${ROOT_BASIC_LIB_LIST}
			${ART_FRAMEWORK_CORE}
		   	${ART_FRAMEWORK_SERVICES_REGISTRY}
		   	${ART_FRAMEWORK_SERVICES_OPTIONAL}
			art_Persistency_Common canvas_Persistency_Common
			art_Persistency_Provenance canvas_Persistency_Provenance
			art_Utilities canvas_Utilities
        )


install_headers()
install_fhicl()
install_source()

