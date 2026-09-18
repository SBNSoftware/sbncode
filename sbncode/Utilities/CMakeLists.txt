cet_make_library(
  INTERFACE
  HEADERS_TARGET_ONLY
  LIBRARIES
    art::Persistency_Common
    canvas::canvas
  )

install_headers()
install_source()

