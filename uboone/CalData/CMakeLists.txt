art_make( BASENAME_ONLY
          LIBRARY_NAME CalData
          LIB_LIBRARIES 
			lardataobj_RawData
                        lardataobj_RecoBase
	                lardata_Utilities
	                ${FHICLCPP}
			${CETLIB}
			${ROOT_BASIC_LIB_LIST}
	  MODULE_LIBRARIES larcore_Geometry
	                CalData
  	  		NoiseFilterAlgs
	  		larcore_Geometry_Geometry_service
	                lardata_Utilities
			larevt_Filters
			lardataobj_RawData
			larevt_CalibrationDBI_IOVData
			larevt_CalibrationDBI_Providers
			lardataobj_RecoBase
			lardata_ArtDataHelper
			SignalShapingServiceMicroBooNE_service
                        ${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
			${ART_FRAMEWORK_SERVICES_REGISTRY}
			${ART_FRAMEWORK_SERVICES_OPTIONAL}
			${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
			art_Persistency_Common canvas_Persistency_Common
			art_Persistency_Provenance canvas_Persistency_Provenance
			art_Utilities canvas_Utilities
			${MF_MESSAGELOGGER}
			${MF_UTILITIES}
	                ${FHICLCPP}
			${CETLIB}
                        ${ROOT_GEOM}
                        ${ROOT_XMLIO}
                        ${ROOT_GDML}
			${ROOT_FFTW}
			${ROOT_BASIC_LIB_LIST}
        )

add_subdirectory(NoiseFilterAlgs)

install_headers()
install_fhicl()
install_source()
