# basic source code CMakeLists.txt

include_directories( $ENV{IFDHC_FQ_DIR}/inc )

art_make(
	LIBRARY_NAME OverlayRawDataMicroBooNE_module
	EXCLUDE RawDigitOverlayMicroBooNEAna_module.cc
	MODULE_LIBRARIES
			ub_RawDigitMixer
			ub_OpDetWaveformMixer
                        lardataobj_RawData
			nusimdata_SimulationBase
			lardataobj_MCBase
			larevt_CalibrationDBI_IOVData
			larevt_CalibrationDBI_Providers
                        ${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_IO_PRODUCTMIX}
			${ART_FRAMEWORK_IO_ROOT}
			${ART_FRAMEWORK_MODULES}
                        ${ART_FRAMEWORK_PRINCIPAL}
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_FRAMEWORK_SERVICES_SYSTEM_FILECATALOGMETADATA_SERVICE}
                        art_Persistency_Common canvas_Persistency_Common
                        art_Utilities canvas_Utilities
			ifdh
                        ${MF_MESSAGELOGGER}
                        ${BOOST_LIB}
                        ${FHICLCPP}
                        ${CETLIB}
)							   

art_make(
	LIBRARY_NAME RawDigitOverlayMicroBooNEAna_module
	EXCLUDE OverlayRawDataMicroBooNE_module.cc
	MODULE_LIBRARIES
			ub_RawDigitAdderAna
                        lardataobj_RawData
                        ${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_IO_ROOT}
			${ART_FRAMEWORK_MODULES}
                        ${ART_FRAMEWORK_PRINCIPAL}
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
			${BOOST_LIB}
                        ${FHICLCPP}
)							   

install_headers()
install_source()
install_fhicl()

# Libraries I didn't need
#			ub_RawDigitAdderAna
#
#			${ART_FRAMEWORK_IO_SOURCES}
#                        art_Persistency_Provenance canvas_Persistency_Provenance
#                        ${MF_UTILITIES}
#                        ${ROOT_BASIC_LIB_LIST}
