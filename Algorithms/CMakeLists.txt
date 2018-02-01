art_make(
          LIB_LIBRARIES
			larcore_Geometry	
		   	uboone_Geometry		
			larcore_Geometry_Geometry_service	
                        ${MF_MESSAGELOGGER}
                        ${MF_UTILITIES}
                        ${CETLIB}
                        ${Boost_SYSTEM_LIBRARY}
                        ${ROOT_BASIC_LIB_LIST}
		   	${ART_FRAMEWORK_SERVICES_REGISTRY}
		   	${ART_FRAMEWORK_SERVICES_OPTIONAL}

	 MODULE_LIBRARIES
			larcore_Geometry
		   	uboone_Geometry
		   	larcore_Geometry_Geometry_service
                        ${ROOT_GEOM}
		   	${ART_FRAMEWORK_SERVICES_REGISTRY}
		   	${ART_FRAMEWORK_SERVICES_OPTIONAL}
        )


install_headers()
install_fhicl()
install_source()

