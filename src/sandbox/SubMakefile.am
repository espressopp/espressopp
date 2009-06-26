if ENABLE_SHARED
python_LTLIBRARIES += sandbox.la
sandbox_la_SOURCES = sandbox/Sandbox.cpp
sandbox_la_LDFLAGS = -shared -module -avoid-version -no-undefined
sandbox_la_LIBADD = $(ESPRESSO_BOOST_LIBS) $(PYTHON_LIBS)
endif
