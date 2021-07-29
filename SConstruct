import os

env = Environment()

env.Append(CCFLAGS='-g',
           CPPPATH=[os.path.join(os.environ['CONDA_PREFIX'],'include'),],
           LIBS=['cantera',],
           LIBPATH=[os.path.join(os.environ['CONDA_PREFIX'],'lib')],
           LINKFLAGS=['-g', '-pthread'])

osap = env.Program('dev-test', ['dev-test.cpp','dev-test-examples.cpp'])

# Default(sample)
