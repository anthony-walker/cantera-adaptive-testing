Flat profile:

Each sample counts as 0.01 seconds.
  %   cumulative   self              self     total
 time   seconds   seconds    calls  ms/call  ms/call  name
 24.76      3.34     3.34   483226     0.01     0.01  void Eigen::SparseMatrix<double, 0, int>::reserveInnerVectors<Eigen::SparseMatrix<double, 0, int>::SingletonVector>(Eigen::SparseMatrix<double, 0, int>::SingletonVector const&)
 13.79      5.20     1.86 2252121019     0.00     0.00  Eigen::internal::CompressedStorage<double, int>::index(long)
 10.90      6.67     1.47 2047437030     0.00     0.00  Eigen::internal::CompressedStorage<double, int>::value(long)
  3.78      7.18     0.51  3715098     0.00     0.00  Eigen::internal::CompressedStorage<double, int>::resize(long, double)
  3.48      7.65     0.47 22299554     0.00     0.00  Eigen::SparseMatrix<double, 0, int>::insertUncompressed(long, long)
  3.34      8.10     0.45                             Cantera::StoichManagerN::dotReaction(double const*)
  3.04      8.51     0.41                             Cantera::AdaptivePreconditioner::solve(Cantera::ReactorNet*, double*, double*)
  2.89      8.90     0.39                             Cantera::AdaptivePreconditioner::TemperatureDerivatives(Cantera::IdealGasConstPressureReactor*, double, double*, double*, double*)
  2.82      9.28     0.38                             Eigen::internal::general_matrix_vector_product<long, double, Eigen::internal::const_blas_data_mapper<double, long, 0>, 0, false, double, Eigen::internal::const_blas_data_mapper<double, long, 0>, false, 0>::run(long, long, Eigen::internal::const_blas_data_mapper<double, long, 0> const&, Eigen::internal::const_blas_data_mapper<double, long, 0> const&, double*, long, double)
  2.08      9.56     0.28                             Cantera::StoichManagerN::multiply(double const*, double*) const
  2.00      9.83     0.27   123948     0.00     0.00  Eigen::internal::CompressedStorage<double, int>::append(double const&, long)
  1.85     10.08     0.25                             Cantera::StoichManagerN::speciesDot(double const*)
  1.48     10.28     0.20                             Eigen::internal::general_matrix_vector_product<long, double, Eigen::internal::const_blas_data_mapper<double, long, 0>, 0, false, double, Eigen::internal::const_blas_data_mapper<double, long, 1>, false, 0>::run(long, long, Eigen::internal::const_blas_data_mapper<double, long, 0> const&, Eigen::internal::const_blas_data_mapper<double, long, 1> const&, double*, long, double)
  1.41     10.47     0.19                             Cantera::NasaPoly1::updateProperties(double const*, double*, double*, double*) const
  1.33     10.65     0.18                             Cantera::MultiBulkRates<Cantera::ArrheniusRate, Cantera::ArrheniusData>::getRateConstants(Cantera::ThermoPhase const&, double*, double*) const
  1.19     10.81     0.16 25528910     0.00     0.00  Eigen::SparseMatrix<double, 0, int>::insert(long, long)
  1.11     10.96     0.15                             void Eigen::internal::sparselu_gemm<double>(long, long, long, double const*, long, double const*, long, double*, long)
  0.82     11.07     0.11                             Eigen::internal::triangular_solve_vector<double, double, long, 1, 5, false, 0>::run(long, double const*, long, double*)
  0.82     11.18     0.11                             Cantera::GasKinetics::updateKc()
  0.82     11.29     0.11                             Cantera::GasKinetics::updateROP()
  0.74     11.39     0.10                             void Eigen::internal::assign_sparse_to_sparse<Eigen::SparseMatrix<double, 0, int>, Eigen::CwiseBinaryOp<Eigen::internal::scalar_difference_op<double, double>, Eigen::SparseMatrix<double, 0, int> const, Eigen::CwiseBinaryOp<Eigen::internal::scalar_product_op<double, double>, Eigen::CwiseNullaryOp<Eigen::internal::scalar_constant_op<double>, Eigen::Matrix<double, -1, -1, 0, -1, -1> const> const, Eigen::SparseMatrix<double, 0, int> const> const> >(Eigen::SparseMatrix<double, 0, int>&, Eigen::CwiseBinaryOp<Eigen::internal::scalar_difference_op<double, double>, Eigen::SparseMatrix<double, 0, int> const, Eigen::CwiseBinaryOp<Eigen::internal::scalar_product_op<double, double>, Eigen::CwiseNullaryOp<Eigen::internal::scalar_constant_op<double>, Eigen::Matrix<double, -1, -1, 0, -1, -1> const> const, Eigen::SparseMatrix<double, 0, int> const> const> const&)
  0.67     11.48     0.09 98094878     0.00     0.00  Eigen::SparseMatrix<double, 0, int>::SingletonVector::operator[](long) const
  0.63     11.57     0.09 98705182     0.00     0.00  int const& std::max<int>(int const&, int const&)
  0.52     11.64     0.07                             int Eigen::internal::coletree<Eigen::SparseMatrix<double, 0, int>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::SparseMatrix<double, 0, int> const&, Eigen::Matrix<int, -1, 1, 0, -1, 1>&, Eigen::Matrix<int, -1, 1, 0, -1, 1>&, Eigen::SparseMatrix<double, 0, int>::StorageIndex*)
  0.52     11.71     0.07 48566472     0.00     0.00  Eigen::SparseCompressedBase<Eigen::SparseMatrix<double, 0, int> >::innerNonZeroPtr() const
  0.44     11.77     0.06                             void Eigen::COLAMDOrdering<int>::operator()<Eigen::SparseMatrix<double, 0, int> >(Eigen::SparseMatrix<double, 0, int> const&, Eigen::PermutationMatrix<-1, -1, int>&)
  0.44     11.83     0.06                             Cantera::GasTransport::getBinDiffCorrection(double, Cantera::MMCollisionInt&, unsigned long, unsigned long, double, double, double&, double&)
  0.44     11.89     0.06                             Cantera::MMCollisionInt::astar(double, double)
  0.37     11.94     0.05                             Eigen::internal::triangular_solve_vector<double, double, long, 1, 2, false, 0>::run(long, double const*, long, double*)
  0.37     11.99     0.05                             void Eigen::internal::conservative_sparse_sparse_product_impl<Eigen::SparseMatrix<double, 0, int>, Eigen::SparseMatrix<double, 0, int>, Eigen::SparseMatrix<double, 0, int> >(Eigen::SparseMatrix<double, 0, int> const&, Eigen::SparseMatrix<double, 0, int> const&, Eigen::SparseMatrix<double, 0, int>&, bool)
  0.37     12.04     0.05                             Cantera::GasKinetics::update_rates_C()
  0.37     12.09     0.05                             Cantera::IdealGasConstPressureReactor::evalEqs(double, double*, double*, double*)
  0.37     12.14     0.05                             N_VDotProd_Serial
  0.37     12.19     0.05                             Cantera::MultiSpeciesThermo::update(double, double*, double*, double*) const
  0.33     12.23     0.05 48566472     0.00     0.00  Eigen::SparseMatrix<double, 0, int>::innerNonZeroPtr() const
  0.30     12.27     0.04                             N_VLinearSum_Serial
  0.30     12.31     0.04                             Cantera::GasKinetics::processFalloffReactions()
  0.30     12.35     0.04                             Cantera::Phase::setMassFractions_NoNorm(double const*)
  0.30     12.39     0.04                             Cantera::Troe::F(double, double const*) const
  0.22     12.42     0.03 26495988     0.00     0.00  Eigen::SparseMatrixBase<Eigen::SparseMatrix<double, 0, int> >::convert_index(long)
  0.22     12.45     0.03 25529536     0.00     0.00  Eigen::SparseMatrix<double, 0, int>::cols() const
  0.22     12.48     0.03   491365     0.00     0.00  Eigen::internal::smart_copy_helper<int, true>::run(int const*, int const*, int*)
  0.22     12.51     0.03                             Eigen::internal::SparseLUImpl<double, int>::column_dfs(long, long, Eigen::Matrix<int, -1, 1, 0, -1, 1>&, long, long&, Eigen::Ref<Eigen::Matrix<int, -1, 1, 0, -1, 1>, 0, Eigen::InnerStride<1> >, Eigen::Matrix<int, -1, 1, 0, -1, 1>&, Eigen::Ref<Eigen::Matrix<int, -1, 1, 0, -1, 1>, 0, Eigen::InnerStride<1> >, Eigen::Matrix<int, -1, 1, 0, -1, 1>&, Eigen::Matrix<int, -1, 1, 0, -1, 1>&, Eigen::Matrix<int, -1, 1, 0, -1, 1>&, Eigen::Matrix<int, -1, 1, 0, -1, 1>&, Eigen::internal::LU_GlobalLU_t<Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<double, -1, 1, 0, -1, 1> >&)
  0.22     12.54     0.03                             Eigen::internal::SparseLUImpl<double, int>::panel_bmod(long, long, long, long, Eigen::Matrix<double, -1, 1, 0, -1, 1>&, Eigen::Matrix<double, -1, 1, 0, -1, 1>&, Eigen::Matrix<int, -1, 1, 0, -1, 1>&, Eigen::Matrix<int, -1, 1, 0, -1, 1>&, Eigen::internal::LU_GlobalLU_t<Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<double, -1, 1, 0, -1, 1> >&)
  0.22     12.57     0.03                             Eigen::internal::SparseLUImpl<double, int>::pivotL(long, double const&, Eigen::Matrix<int, -1, 1, 0, -1, 1>&, Eigen::Matrix<int, -1, 1, 0, -1, 1>&, long&, Eigen::internal::LU_GlobalLU_t<Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<double, -1, 1, 0, -1, 1> >&)
  0.22     12.60     0.03                             Eigen::internal::SparseLUImpl<double, int>::panel_dfs(long, long, long, Eigen::SparseMatrix<double, 0, int>&, Eigen::Matrix<int, -1, 1, 0, -1, 1>&, long&, Eigen::Matrix<double, -1, 1, 0, -1, 1>&, Eigen::Matrix<int, -1, 1, 0, -1, 1>&, Eigen::Matrix<int, -1, 1, 0, -1, 1>&, Eigen::Matrix<int, -1, 1, 0, -1, 1>&, Eigen::Matrix<int, -1, 1, 0, -1, 1>&, Eigen::Matrix<int, -1, 1, 0, -1, 1>&, Eigen::Matrix<int, -1, 1, 0, -1, 1>&, Eigen::Matrix<int, -1, 1, 0, -1, 1>&, Eigen::internal::LU_GlobalLU_t<Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<double, -1, 1, 0, -1, 1> >&)
  0.22     12.63     0.03                             Cantera::AdaptivePreconditioner::getPreconditionerType()
  0.22     12.66     0.03                             Cantera::polyfit(unsigned long, unsigned long, double const*, double const*, double const*, double*)
  0.22     12.69     0.03                             Cantera::Kinetics::getRevReactionDelta(double const*, double*)
  0.22     12.72     0.03                             void Eigen::internal::MappedSuperNodalMatrix<double, int>::solveInPlace<Eigen::Matrix<double, -1, 1, 0, -1, 1> >(Eigen::MatrixBase<Eigen::Matrix<double, -1, 1, 0, -1, 1> >&) const
  0.22     12.75     0.03                             Cantera::Phase::mean_X(std::vector<double, std::allocator<double> > const&) const
  0.15     12.77     0.02 51059698     0.00     0.00  Eigen::internal::CompressedStorage<double, int>::allocatedSize() const
  0.15     12.79     0.02 48565846     0.00     0.00  Eigen::SparseCompressedBase<Eigen::SparseMatrix<double, 0, int> >::isCompressed() const
  0.15     12.81     0.02 25529536     0.00     0.00  Eigen::SparseMatrix<double, 0, int>::rows() const
  0.15     12.83     0.02   491365     0.00     0.00  void Eigen::internal::smart_copy<double>(double const*, double const*, double*)
  0.15     12.85     0.02      626     0.03     0.08  void Eigen::SparseMatrix<double, 0, int>::reserveInnerVectors<Eigen::CwiseNullaryOp<Eigen::internal::scalar_constant_op<int>, Eigen::Array<int, -1, 1, 0, -1, 1> > >(Eigen::CwiseNullaryOp<Eigen::internal::scalar_constant_op<int>, Eigen::Array<int, -1, 1, 0, -1, 1> > const&)
  0.15     12.87     0.02                             N_VWSqrSumLocal_Serial
  0.15     12.89     0.02                             Eigen::internal::binary_evaluator<Eigen::CwiseBinaryOp<Eigen::internal::scalar_difference_op<double, double>, Eigen::SparseMatrix<double, 0, int> const, Eigen::SparseMatrix<double, 0, int> const>, Eigen::internal::IteratorBased, Eigen::internal::IteratorBased, double, double>::InnerIterator::operator++()
  0.15     12.91     0.02                             void Eigen::internal::assign_sparse_to_sparse<Eigen::SparseMatrix<double, 0, int>, Eigen::Map<Eigen::SparseMatrix<double, 0, int>, 0, Eigen::Stride<0, 0> > >(Eigen::SparseMatrix<double, 0, int>&, Eigen::Map<Eigen::SparseMatrix<double, 0, int>, 0, Eigen::Stride<0, 0> > const&)
  0.15     12.93     0.02                             Cantera::quadInterp(double, double*, double*)
  0.15     12.95     0.02                             Cantera::GasKinetics::update_rates_T()
  0.15     12.97     0.02                             Cantera::GasTransport::fitDiffCoeffs(Cantera::MMCollisionInt&)
  0.15     12.99     0.02                             Cantera::IdealGasPhase::getPartialMolarEnthalpies(double*) const
  0.15     13.01     0.02                             Cantera::Troe::updateTemp(double, double*) const
  0.15     13.03     0.02                             Cantera::Phase::saveState(unsigned long, double*) const
  0.15     13.05     0.02                             Cantera::NasaPoly2::temperaturePolySize() const
  0.11     13.07     0.02 26619937     0.00     0.00  int Eigen::internal::convert_index<int, long>(long const&)
  0.07     13.08     0.01  7101344     0.00     0.00  Eigen::SparseCompressedBase<Eigen::SparseMatrix<double, 0, int> >::InnerIterator::operator bool() const
  0.07     13.09     0.01  6849692     0.00     0.00  Eigen::SparseCompressedBase<Eigen::SparseMatrix<double, 0, int> >::InnerIterator::index() const
  0.07     13.10     0.01  3555680     0.00     0.00  Eigen::internal::variable_if_dynamic<long, 0>::variable_if_dynamic(long)
  0.07     13.11     0.01  3550672     0.00     0.00  Eigen::DenseCoeffsBase<Eigen::Matrix<int, -1, 1, 0, -1, 1>, 1>::operator[](long)
  0.07     13.12     0.01   982124     0.00     0.00  std::enable_if<std::__and_<std::__not_<std::__is_tuple_like<int*> >, std::is_move_constructible<int*>, std::is_move_assignable<int*> >::value, void>::type std::swap<int*>(int*&, int*&)
  0.07     13.13     0.01   504556     0.00     0.00  Eigen::SparseCompressedBase<Eigen::SparseMatrix<double, 0, int> >::outerIndexPtr() const
  0.07     13.14     0.01   127078     0.00     0.00  Eigen::internal::evaluator_base<Eigen::CwiseNullaryOp<Eigen::internal::scalar_constant_op<int>, Eigen::Array<int, -1, 1, 0, -1, 1> > >::~evaluator_base()
  0.07     13.15     0.01   127078     0.00     0.00  int Eigen::internal::evaluator<Eigen::CwiseNullaryOp<Eigen::internal::scalar_constant_op<int>, Eigen::Array<int, -1, 1, 0, -1, 1> > >::coeff<long>(long) const
  0.07     13.16     0.01    31300     0.00     0.00  long long __vector(2) Eigen::internal::pset1<long long __vector(2)>(Eigen::internal::unpacket_traits<long long __vector(2)>::type const&)
  0.07     13.17     0.01     1253     0.01     0.01  Eigen::SparseMatrix<double, 0, int>::SparseMatrix(long, long)
  0.07     13.18     0.01     1252     0.01     0.01  Eigen::SparseCompressedBase<Eigen::SparseMatrix<double, 0, int> >::nonZeros() const
  0.07     13.19     0.01                             N_VProd_Serial
  0.07     13.20     0.01                             N_VSetArrayPointer_Serial
  0.07     13.21     0.01                             N_VWrmsNorm_Serial
  0.07     13.22     0.01                             YAML::NodeBuilder::OnMapStart(YAML::Mark const&, std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > const&, unsigned long, YAML::EmitterStyle::value)
  0.07     13.23     0.01                             Eigen::internal::SparseLUImpl<double, int>::column_bmod(long, long, Eigen::Ref<Eigen::Matrix<double, -1, 1, 0, -1, 1>, 0, Eigen::InnerStride<1> >, Eigen::Matrix<double, -1, 1, 0, -1, 1>&, Eigen::Ref<Eigen::Matrix<int, -1, 1, 0, -1, 1>, 0, Eigen::InnerStride<1> >, Eigen::Ref<Eigen::Matrix<int, -1, 1, 0, -1, 1>, 0, Eigen::InnerStride<1> >, long, Eigen::internal::LU_GlobalLU_t<Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<double, -1, 1, 0, -1, 1> >&)
  0.07     13.24     0.01                             Eigen::internal::SparseLUImpl<double, int>::pruneL(long, Eigen::Matrix<int, -1, 1, 0, -1, 1> const&, long, long, Eigen::Matrix<int, -1, 1, 0, -1, 1> const&, Eigen::Ref<Eigen::Matrix<int, -1, 1, 0, -1, 1>, 0, Eigen::InnerStride<1> >, Eigen::Matrix<int, -1, 1, 0, -1, 1>&, Eigen::internal::LU_GlobalLU_t<Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<double, -1, 1, 0, -1, 1> >&)
  0.07     13.25     0.01                             void Eigen::internal::LU_kernel_bmod<1>::run<Eigen::VectorBlock<Eigen::Matrix<double, -1, 1, 0, -1, 1>, -1>, Eigen::Matrix<double, -1, 1, 0, -1, 1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(long, Eigen::VectorBlock<Eigen::Matrix<double, -1, 1, 0, -1, 1>, -1>&, Eigen::Matrix<double, -1, 1, 0, -1, 1>&, Eigen::Matrix<double, -1, 1, 0, -1, 1>&, long&, long, long, Eigen::Matrix<int, -1, 1, 0, -1, 1>&, long, long)
  0.07     13.26     0.01                             void Eigen::internal::assign_sparse_to_sparse<Eigen::SparseMatrix<double, 0, int>, Eigen::Product<Eigen::SparseMatrix<double, 0, int>, Eigen::CwiseBinaryOp<Eigen::internal::scalar_difference_op<double, double>, Eigen::SparseMatrix<double, 0, int> const, Eigen::SparseMatrix<double, 0, int> const>, 2> >(Eigen::SparseMatrix<double, 0, int>&, Eigen::Product<Eigen::SparseMatrix<double, 0, int>, Eigen::CwiseBinaryOp<Eigen::internal::scalar_difference_op<double, double>, Eigen::SparseMatrix<double, 0, int> const, Eigen::SparseMatrix<double, 0, int> const>, 2> const&)
  0.07     13.27     0.01                             Eigen::internal::triangular_solve_matrix<double, long, 1, 5, false, 0, 0>::run(long, long, double const*, long, double*, long, Eigen::internal::level3_blocking<double, double>&)
  0.07     13.28     0.01                             Cantera::FlowDevice::outletSpeciesMassFlowRate(unsigned long)
  0.07     13.29     0.01                             Cantera::MMCollisionInt::bstar(double, double)
  0.07     13.30     0.01                             Cantera::MMCollisionInt::cstar(double, double)
  0.07     13.31     0.01                             Cantera::MMCollisionInt::omega22(double, double)
  0.07     13.32     0.01                             Cantera::AdaptivePreconditioner::ForwardConversion(Cantera::Reactor*, double*, double*, unsigned long)
  0.07     13.33     0.01                             Cantera::AdaptivePreconditioner::preconditionerErrorCheck()
  0.07     13.34     0.01                             Cantera::Reactor::updateConnected(bool)
  0.07     13.35     0.01                             Cantera::Kinetics::getNetProductionRates(double*)
  0.07     13.36     0.01                             Cantera::IdealGasPhase::getPureGibbs(double*) const
  0.07     13.37     0.01                             Cantera::IdealGasPhase::pressure() const
  0.07     13.38     0.01                             Cantera::ThreeBodyReaction3::type[abi:cxx11]() const
  0.07     13.39     0.01                             Cantera::SRI::getParameters(double*) const
  0.07     13.40     0.01                             Cantera::Kinetics::reactionTypeStr[abi:cxx11](unsigned long) const
  0.07     13.41     0.01                             Cantera::NasaPoly2::updateProperties(double const*, double*, double*, double*) const
  0.07     13.42     0.01                             std::map<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, unsigned long, std::less<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > >, std::allocator<std::pair<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > const, unsigned long> > >::at(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > const&)
  0.07     13.43     0.01                             std::vector<Cantera::AnyMap, std::allocator<Cantera::AnyMap> >::vector(std::vector<Cantera::AnyMap, std::allocator<Cantera::AnyMap> > const&)
  0.07     13.44     0.01                             void std::vector<std::vector<double, std::allocator<double> >, std::allocator<std::vector<double, std::allocator<double> > > >::_M_realloc_insert<std::vector<double, std::allocator<double> > >(__gnu_cxx::__normal_iterator<std::vector<double, std::allocator<double> >*, std::vector<std::vector<double, std::allocator<double> >, std::allocator<std::vector<double, std::allocator<double> > > > >, std::vector<double, std::allocator<double> >&&)
  0.07     13.45     0.01                             std::_Rb_tree<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, std::pair<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > const, unsigned long>, std::_Select1st<std::pair<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > const, unsigned long> >, std::less<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > >, std::allocator<std::pair<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > const, unsigned long> > >::_M_get_insert_unique_pos(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > const&)
  0.07     13.46     0.01                             cvLsSolve
  0.07     13.47     0.01                             cvStep
  0.07     13.48     0.01                             cvodes_err
  0.04     13.48     0.01   496394     0.00     0.00  Eigen::internal::CompressedStorage<double, int>::reallocate(long)
  0.04     13.49     0.01        2     2.50     2.50  std::__shared_count<(__gnu_cxx::_Lock_policy)2>::__shared_count(std::__shared_count<(__gnu_cxx::_Lock_policy)2> const&)
  0.04     13.49     0.01                             unsigned long* std::uninitialized_copy<unsigned long const*, unsigned long*>(unsigned long const*, unsigned long const*, unsigned long*)
  0.00     13.49     0.00 49578088     0.00     0.00  Eigen::SparseMatrixBase<Eigen::SparseMatrix<double, 0, int> >::derived() const
  0.00     13.49     0.00 28752011     0.00     0.00  Eigen::internal::CompressedStorage<double, int>::size() const
  0.00     13.49     0.00  7101344     0.00     0.00  Eigen::EigenBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >::derived() const
  0.00     13.49     0.00  6849692     0.00     0.00  Eigen::SparseCompressedBase<Eigen::SparseMatrix<double, 0, int> >::InnerIterator::operator++()
  0.00     13.49     0.00  4673043     0.00     0.00  Eigen::internal::noncopyable::noncopyable()
  0.00     13.49     0.00  4673043     0.00     0.00  Eigen::internal::noncopyable::~noncopyable()
  0.00     13.49     0.00  3550672     0.00     0.00  Eigen::DenseStorage<int, -1, -1, 1, 0>::cols()
  0.00     13.49     0.00  3550672     0.00     0.00  Eigen::DenseCoeffsBase<Eigen::Matrix<int, -1, 1, 0, -1, 1>, 1>::coeffRef(long)
  0.00     13.49     0.00  3550672     0.00     0.00  Eigen::internal::evaluator_base<Eigen::Matrix<int, -1, 1, 0, -1, 1> >::evaluator_base()
  0.00     13.49     0.00  3550672     0.00     0.00  Eigen::internal::evaluator_base<Eigen::Matrix<int, -1, 1, 0, -1, 1> >::~evaluator_base()
  0.00     13.49     0.00  3550672     0.00     0.00  Eigen::internal::evaluator<Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> > >::coeffRef(long)
  0.00     13.49     0.00  3550672     0.00     0.00  Eigen::internal::evaluator<Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> > >::evaluator(Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> > const&)
  0.00     13.49     0.00  3550672     0.00     0.00  Eigen::internal::evaluator<Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> > >::~evaluator()
  0.00     13.49     0.00  3550672     0.00     0.00  Eigen::internal::evaluator<Eigen::Matrix<int, -1, 1, 0, -1, 1> >::evaluator(Eigen::Matrix<int, -1, 1, 0, -1, 1> const&)
  0.00     13.49     0.00  3550672     0.00     0.00  Eigen::internal::evaluator<Eigen::Matrix<int, -1, 1, 0, -1, 1> >::~evaluator()
  0.00     13.49     0.00  3550672     0.00     0.00  Eigen::EigenBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >::derived()
  0.00     13.49     0.00  3550672     0.00     0.00  Eigen::DenseStorage<int, -1, -1, 1, 0>::data() const
  0.00     13.49     0.00  3550672     0.00     0.00  Eigen::DenseStorage<int, -1, -1, 1, 0>::rows() const
  0.00     13.49     0.00  3550672     0.00     0.00  Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >::cols() const
  0.00     13.49     0.00  3550672     0.00     0.00  Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >::data() const
  0.00     13.49     0.00  3550672     0.00     0.00  Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >::rows() const
  0.00     13.49     0.00  3550672     0.00     0.00  Eigen::EigenBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >::cols() const
  0.00     13.49     0.00  3550672     0.00     0.00  Eigen::EigenBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >::rows() const
  0.00     13.49     0.00  3550672     0.00     0.00  Eigen::EigenBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >::size() const
  0.00     13.49     0.00  3424846     0.00     0.00  Eigen::SparseCompressedBase<Eigen::SparseMatrix<double, 0, int> >::InnerIterator::value() const
  0.00     13.49     0.00  2946372     0.00     0.00  std::remove_reference<int*&>::type&& std::move<int*&>(int*&)
  0.00     13.49     0.00  1491060     0.00     0.00  std::remove_reference<double*&>::type&& std::move<double*&>(double*&)
  0.00     13.49     0.00   987759     0.00     0.00  Eigen::internal::scoped_array<double>::ptr()
  0.00     13.49     0.00   987759     0.00     0.00  Eigen::internal::scoped_array<int>::ptr()
  0.00     13.49     0.00   981504     0.00     0.00  long const& std::min<long>(long const&, long const&)
  0.00     13.49     0.00   504556     0.00     0.00  Eigen::SparseMatrix<double, 0, int>::outerIndexPtr() const
  0.00     13.49     0.00   497020     0.00     0.00  std::enable_if<std::__and_<std::__not_<std::__is_tuple_like<double*> >, std::is_move_constructible<double*>, std::is_move_assignable<double*> >::value, void>::type std::swap<double*>(double*&, double*&)
  0.00     13.49     0.00   496394     0.00     0.00  Eigen::internal::scoped_array<double>::scoped_array(long)
  0.00     13.49     0.00   496394     0.00     0.00  Eigen::internal::scoped_array<double>::~scoped_array()
  0.00     13.49     0.00   496394     0.00     0.00  Eigen::internal::scoped_array<int>::scoped_array(long)
  0.00     13.49     0.00   496394     0.00     0.00  Eigen::internal::scoped_array<int>::~scoped_array()
  0.00     13.49     0.00   491365     0.00     0.00  void Eigen::internal::smart_copy<int>(int const*, int const*, int*)
  0.00     13.49     0.00   491365     0.00     0.00  Eigen::internal::smart_copy_helper<double, true>::run(double const*, double const*, double*)
  0.00     13.49     0.00   485110     0.00     0.00  Eigen::GenericNumTraits<int>::highest()
  0.00     13.49     0.00   485110     0.00     0.00  std::numeric_limits<int>::max()
  0.00     13.49     0.00   483226     0.00     0.00  Eigen::SparseMatrix<double, 0, int>::SingletonVector::SingletonVector(long, long)
  0.00     13.49     0.00   483226     0.00     0.01  void Eigen::SparseMatrix<double, 0, int>::reserve<Eigen::SparseMatrix<double, 0, int>::SingletonVector>(Eigen::SparseMatrix<double, 0, int>::SingletonVector const&, Eigen::SparseMatrix<double, 0, int>::SingletonVector::value_type const&)
  0.00     13.49     0.00   483226     0.00     0.00  void Eigen::internal::ignore_unused_variable<int>(int const&)
  0.00     13.49     0.00   381234     0.00     0.00  Eigen::EigenBase<Eigen::CwiseNullaryOp<Eigen::internal::scalar_constant_op<int>, Eigen::Array<int, -1, 1, 0, -1, 1> > >::derived() const
  0.00     13.49     0.00   255408     0.00     0.00  Eigen::internal::variable_if_dynamic<long, -1>::variable_if_dynamic(long)
  0.00     13.49     0.00   254782     0.00     0.00  Eigen::SparseMatrix<double, 0, int>::outerSize() const
  0.00     13.49     0.00   251652     0.00     0.00  Eigen::SparseCompressedBase<Eigen::SparseMatrix<double, 0, int> >::InnerIterator::InnerIterator(Eigen::SparseCompressedBase<Eigen::SparseMatrix<double, 0, int> > const&, long)
  0.00     13.49     0.00   251652     0.00     0.00  Eigen::internal::evaluator<Eigen::SparseCompressedBase<Eigen::SparseMatrix<double, 0, int> > >::operator Eigen::SparseMatrix<double, 0, int>&()
  0.00     13.49     0.00   251652     0.00     0.00  Eigen::SparseMatrix<double, 0, int>::innerIndexPtr() const
  0.00     13.49     0.00   251652     0.00     0.00  Eigen::SparseMatrix<double, 0, int>::valuePtr() const
  0.00     13.49     0.00   251652     0.00     0.00  Eigen::SparseMatrixBase<Eigen::SparseMatrix<double, 0, int> >::const_cast_derived() const
  0.00     13.49     0.00   251652     0.00     0.00  Eigen::SparseCompressedBase<Eigen::SparseMatrix<double, 0, int> >::innerIndexPtr() const
  0.00     13.49     0.00   251652     0.00     0.00  Eigen::SparseCompressedBase<Eigen::SparseMatrix<double, 0, int> >::valuePtr() const
  0.00     13.49     0.00   251652     0.00     0.00  Eigen::internal::CompressedStorage<double, int>::indexPtr() const
  0.00     13.49     0.00   251652     0.00     0.00  Eigen::internal::CompressedStorage<double, int>::valuePtr() const
  0.00     13.49     0.00   232145     0.00     0.00  __gthread_active_p()
  0.00     13.49     0.00   232143     0.00     0.00  __gnu_cxx::__exchange_and_add(int volatile*, int)
  0.00     13.49     0.00   232143     0.00     0.00  __gnu_cxx::__exchange_and_add_dispatch(int*, int)
  0.00     13.49     0.00   221991     0.00     0.00  std::_Sp_counted_base<(__gnu_cxx::_Lock_policy)2>::_M_release()
  0.00     13.49     0.00   199068     0.00     0.00  Eigen::internal::variable_if_dynamic<long, 1>::value()
  0.00     13.49     0.00   133964     0.00     0.00  Eigen::internal::variable_if_dynamic<long, -1>::value() const
  0.00     13.49     0.00   128956     0.00     0.00  Eigen::internal::scalar_constant_op<int>::scalar_constant_op(Eigen::internal::scalar_constant_op<int> const&)
  0.00     13.49     0.00   128330     0.00     0.00  Eigen::SparseMatrix<double, 1, int>::outerSize() const
  0.00     13.49     0.00   127704     0.00     0.00  int Eigen::internal::nullary_wrapper<int, Eigen::internal::scalar_constant_op<int>, true, false, false>::operator()<long>(Eigen::internal::scalar_constant_op<int> const&, long, long) const
  0.00     13.49     0.00   127704     0.00     0.00  Eigen::internal::scalar_constant_op<int>::operator()() const
  0.00     13.49     0.00   127078     0.00     0.00  Eigen::internal::evaluator_base<Eigen::CwiseNullaryOp<Eigen::internal::scalar_constant_op<int>, Eigen::Array<int, -1, 1, 0, -1, 1> > >::evaluator_base()
  0.00     13.49     0.00   127078     0.00     0.00  Eigen::internal::evaluator<Eigen::CwiseNullaryOp<Eigen::internal::scalar_constant_op<int>, Eigen::Array<int, -1, 1, 0, -1, 1> > >::evaluator(Eigen::CwiseNullaryOp<Eigen::internal::scalar_constant_op<int>, Eigen::Array<int, -1, 1, 0, -1, 1> > const&)
  0.00     13.49     0.00   127078     0.00     0.00  Eigen::internal::evaluator<Eigen::CwiseNullaryOp<Eigen::internal::scalar_constant_op<int>, Eigen::Array<int, -1, 1, 0, -1, 1> > >::~evaluator()
  0.00     13.49     0.00   127078     0.00     0.00  Eigen::CwiseNullaryOp<Eigen::internal::scalar_constant_op<int>, Eigen::Array<int, -1, 1, 0, -1, 1> >::cols() const
  0.00     13.49     0.00   127078     0.00     0.00  Eigen::CwiseNullaryOp<Eigen::internal::scalar_constant_op<int>, Eigen::Array<int, -1, 1, 0, -1, 1> >::rows() const
  0.00     13.49     0.00   127078     0.00     0.00  Eigen::CwiseNullaryOp<Eigen::internal::scalar_constant_op<int>, Eigen::Array<int, -1, 1, 0, -1, 1> >::functor() const
  0.00     13.49     0.00   127078     0.00     0.00  Eigen::DenseCoeffsBase<Eigen::CwiseNullaryOp<Eigen::internal::scalar_constant_op<int>, Eigen::Array<int, -1, 1, 0, -1, 1> >, 0>::coeff(long) const
  0.00     13.49     0.00   127078     0.00     0.00  Eigen::DenseCoeffsBase<Eigen::CwiseNullaryOp<Eigen::internal::scalar_constant_op<int>, Eigen::Array<int, -1, 1, 0, -1, 1> >, 0>::operator[](long) const
  0.00     13.49     0.00   127078     0.00     0.00  Eigen::EigenBase<Eigen::CwiseNullaryOp<Eigen::internal::scalar_constant_op<int>, Eigen::Array<int, -1, 1, 0, -1, 1> > >::cols() const
  0.00     13.49     0.00   127078     0.00     0.00  Eigen::EigenBase<Eigen::CwiseNullaryOp<Eigen::internal::scalar_constant_op<int>, Eigen::Array<int, -1, 1, 0, -1, 1> > >::rows() const
  0.00     13.49     0.00   127078     0.00     0.00  Eigen::EigenBase<Eigen::CwiseNullaryOp<Eigen::internal::scalar_constant_op<int>, Eigen::Array<int, -1, 1, 0, -1, 1> > >::size() const
  0.00     13.49     0.00    74974     0.00     0.00  Eigen::internal::aligned_malloc(unsigned long)
  0.00     13.49     0.00    74974     0.00     0.00  Eigen::internal::check_that_malloc_is_allowed()
  0.00     13.49     0.00    31926     0.00     0.00  Eigen::internal::mapbase_evaluator<Eigen::Map<Eigen::Matrix<int, -1, 1, 0, -1, 1>, 0, Eigen::Stride<0, 0> >, Eigen::Matrix<int, -1, 1, 0, -1, 1> >::coeffRef(long)
  0.00     13.49     0.00    31300     0.00     0.00  void Eigen::internal::generic_dense_assignment_kernel<Eigen::internal::evaluator<Eigen::Map<Eigen::Matrix<int, -1, 1, 0, -1, 1>, 0, Eigen::Stride<0, 0> > >, Eigen::internal::evaluator<Eigen::CwiseNullaryOp<Eigen::internal::scalar_constant_op<int>, Eigen::Matrix<int, -1, 1, 0, -1, 1> > >, Eigen::internal::assign_op<int, int>, 0>::assignPacket<16, 0, long long __vector(2)>(long)
  0.00     13.49     0.00    31300     0.00     0.00  long long __vector(2) Eigen::internal::ploadu<long long __vector(2)>(Eigen::internal::unpacket_traits<long long __vector(2)>::type const*)
  0.00     13.49     0.00    31300     0.00     0.00  void Eigen::internal::pstore<int, long long __vector(2)>(int*, long long __vector(2) const&)
  0.00     13.49     0.00    31300     0.00     0.00  long long __vector(2) Eigen::internal::nullary_wrapper<int, Eigen::internal::scalar_constant_op<int>, true, false, false>::packetOp<long long __vector(2), long>(Eigen::internal::scalar_constant_op<int> const&, long, long) const
  0.00     13.49     0.00    31300     0.00     0.00  long long __vector(2) Eigen::internal::redux_evaluator<Eigen::Map<Eigen::Matrix<int, -1, 1, 0, -1, 1> const, 0, Eigen::Stride<0, 0> > >::packet<0, long long __vector(2)>(long) const
  0.00     13.49     0.00    31300     0.00     0.00  long long __vector(2) Eigen::internal::mapbase_evaluator<Eigen::Map<Eigen::Matrix<int, -1, 1, 0, -1, 1> const, 0, Eigen::Stride<0, 0> >, Eigen::Matrix<int, -1, 1, 0, -1, 1> const>::packet<0, long long __vector(2)>(long) const
  0.00     13.49     0.00    31300     0.00     0.00  long long __vector(2) const Eigen::internal::scalar_constant_op<int>::packetOp<long long __vector(2)>() const
  0.00     13.49     0.00    31300     0.00     0.00  void Eigen::internal::assign_op<int, int>::assignPacket<16, long long __vector(2)>(int*, long long __vector(2) const&) const
  0.00     13.49     0.00    31300     0.00     0.00  long long __vector(2) Eigen::internal::evaluator<Eigen::CwiseNullaryOp<Eigen::internal::scalar_constant_op<int>, Eigen::Matrix<int, -1, 1, 0, -1, 1> > >::packet<0, long long __vector(2), long>(long) const
  0.00     13.49     0.00    30674     0.00     0.00  long long __vector(2) Eigen::internal::padd<long long __vector(2)>(long long __vector(2) const&, long long __vector(2) const&)
  0.00     13.49     0.00    30674     0.00     0.00  long long __vector(2) const Eigen::internal::scalar_sum_op<int, int>::packetOp<long long __vector(2)>(long long __vector(2) const&, long long __vector(2) const&) const
  0.00     13.49     0.00     9400     0.00     0.00  Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >::_check_template_params()
  0.00     13.49     0.00     8776     0.00     0.00  Eigen::internal::aligned_free(void*)
  0.00     13.49     0.00     8776     0.00     0.00  void Eigen::internal::conditional_aligned_free<true>(void*)
  0.00     13.49     0.00     8774     0.00     0.00  void Eigen::internal::conditional_aligned_delete_auto<int, true>(int*, unsigned long)
  0.00     13.49     0.00     8148     0.00     0.00  Eigen::MatrixBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >::MatrixBase()
  0.00     13.49     0.00     8148     0.00     0.00  Eigen::DenseStorage<int, -1, -1, 1, 0>::DenseStorage()
  0.00     13.49     0.00     8148     0.00     0.00  Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >::PlainObjectBase()
  0.00     13.49     0.00     8148     0.00     0.00  Eigen::DenseBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >::DenseBase()
  0.00     13.49     0.00     8138     0.00     0.00  Eigen::DenseStorage<int, -1, -1, 1, 0>::resize(long, long, long)
  0.00     13.49     0.00     8138     0.00     0.00  void Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >::_init1<long>(long, Eigen::internal::enable_if<((((Eigen::DenseBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >::{unnamed type#1})-1)!=(1))||(!Eigen::internal::is_convertible<long, int>::value))&&((!((Eigen::internal::is_same<Eigen::MatrixXpr, Eigen::ArrayXpr>::{unnamed type#1})0))||((({unnamed type#1})-1)==Eigen::Dynamic)), Eigen::internal::is_convertible>::type*)
  0.00     13.49     0.00     8138     0.00     0.00  Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >::resize(long)
  0.00     13.49     0.00     8138     0.00     0.00  Eigen::Matrix<int, -1, 1, 0, -1, 1>::Matrix<long>(long const&)
  0.00     13.49     0.00     8138     0.00     0.00  void Eigen::internal::ignore_unused_variable<bool>(bool const&)
  0.00     13.49     0.00     8138     0.00     0.00  void* Eigen::internal::conditional_aligned_malloc<true>(unsigned long)
  0.00     13.49     0.00     8138     0.00     0.00  int* Eigen::internal::conditional_aligned_new_auto<int, true>(unsigned long)
  0.00     13.49     0.00     8138     0.00     0.00  Eigen::EigenBase<Eigen::Map<Eigen::Matrix<int, -1, 1, 0, -1, 1> const, 0, Eigen::Stride<0, 0> > >::derived() const
  0.00     13.49     0.00     7512     0.00     0.00  std::remove_reference<long&>::type&& std::move<long&>(long&)
  0.00     13.49     0.00     6904     0.00     0.00  Eigen::internal::CompressedStorage<double, int>::clear()
  0.00     13.49     0.00     6278     0.00     0.00  Eigen::SparseMatrix<double, 0, int>::resize(long, long)
  0.00     13.49     0.00     5651     0.00     0.00  Eigen::internal::CompressedStorage<double, int>::~CompressedStorage()
  0.00     13.49     0.00     5025     0.00     0.00  Eigen::SparseMatrix<double, 0, int>::~SparseMatrix()
  0.00     13.49     0.00     3756     0.00     0.00  Eigen::internal::variable_if_dynamic<long, 1>::variable_if_dynamic(long)
  0.00     13.49     0.00     3756     0.00     0.00  Eigen::EigenBase<Eigen::Map<Eigen::Matrix<int, -1, 1, 0, -1, 1>, 0, Eigen::Stride<0, 0> > >::derived() const
  0.00     13.49     0.00     3130     0.00     0.00  Eigen::MapBase<Eigen::Map<Eigen::Matrix<int, -1, 1, 0, -1, 1> const, 0, Eigen::Stride<0, 0> >, 0>::cols() const
  0.00     13.49     0.00     3130     0.00     0.00  Eigen::MapBase<Eigen::Map<Eigen::Matrix<int, -1, 1, 0, -1, 1> const, 0, Eigen::Stride<0, 0> >, 0>::rows() const
  0.00     13.49     0.00     3130     0.00     0.00  Eigen::MapBase<Eigen::Map<Eigen::Matrix<int, -1, 1, 0, -1, 1>, 0, Eigen::Stride<0, 0> >, 0>::cols() const
  0.00     13.49     0.00     3130     0.00     0.00  Eigen::MapBase<Eigen::Map<Eigen::Matrix<int, -1, 1, 0, -1, 1>, 0, Eigen::Stride<0, 0> >, 0>::rows() const
  0.00     13.49     0.00     3130     0.00     0.00  Eigen::EigenBase<Eigen::Map<Eigen::Matrix<int, -1, 1, 0, -1, 1> const, 0, Eigen::Stride<0, 0> > >::cols() const
  0.00     13.49     0.00     3130     0.00     0.00  Eigen::EigenBase<Eigen::Map<Eigen::Matrix<int, -1, 1, 0, -1, 1> const, 0, Eigen::Stride<0, 0> > >::rows() const
  0.00     13.49     0.00     2504     0.00     0.00  Eigen::internal::variable_if_dynamic<long, 0>::value()
  0.00     13.49     0.00     2504     0.00     0.00  Eigen::EigenBase<Eigen::Map<Eigen::Matrix<int, -1, 1, 0, -1, 1> const, 0, Eigen::Stride<0, 0> > >::size() const
  0.00     13.49     0.00     2504     0.00     0.00  std::enable_if<std::__and_<std::__not_<std::__is_tuple_like<long> >, std::is_move_constructible<long>, std::is_move_assignable<long> >::value, void>::type std::swap<long>(long&, long&)
  0.00     13.49     0.00     1882     0.00     0.00  Eigen::internal::CompressedStorage<double, int>::CompressedStorage()
  0.00     13.49     0.00     1878     0.00     0.00  Eigen::EigenBase<Eigen::Map<Eigen::Matrix<int, -1, 1, 0, -1, 1>, 0, Eigen::Stride<0, 0> > >::derived()
  0.00     13.49     0.00     1878     0.00     0.00  Eigen::internal::scalar_sum_op<int, int>::operator()(int const&, int const&) const
  0.00     13.49     0.00     1878     0.00     0.00  Eigen::internal::redux_evaluator<Eigen::Map<Eigen::Matrix<int, -1, 1, 0, -1, 1> const, 0, Eigen::Stride<0, 0> > >::coeff(long) const
  0.00     13.49     0.00     1878     0.00     0.00  Eigen::internal::mapbase_evaluator<Eigen::Map<Eigen::Matrix<int, -1, 1, 0, -1, 1> const, 0, Eigen::Stride<0, 0> >, Eigen::Matrix<int, -1, 1, 0, -1, 1> const>::coeff(long) const
  0.00     13.49     0.00     1878     0.00     0.00  Eigen::EigenBase<Eigen::Map<Eigen::Matrix<int, -1, 1, 0, -1, 1>, 0, Eigen::Stride<0, 0> > >::cols() const
  0.00     13.49     0.00     1878     0.00     0.00  Eigen::EigenBase<Eigen::Map<Eigen::Matrix<int, -1, 1, 0, -1, 1>, 0, Eigen::Stride<0, 0> > >::rows() const
  0.00     13.49     0.00     1256     0.00     0.00  Eigen::SparseMatrix<double, 0, int>::check_template_parameters()
  0.00     13.49     0.00     1256     0.00     0.00  Eigen::SparseMatrixBase<Eigen::SparseMatrix<double, 0, int> >::SparseMatrixBase()
  0.00     13.49     0.00     1256     0.00     0.00  Eigen::SparseCompressedBase<Eigen::SparseMatrix<double, 0, int> >::SparseCompressedBase()
  0.00     13.49     0.00     1252     0.00     0.00  Eigen::Stride<0, 0>::Stride(Eigen::Stride<0, 0> const&)
  0.00     13.49     0.00     1252     0.00     0.00  Eigen::Stride<0, 0>::Stride()
  0.00     13.49     0.00     1252     0.00     0.00  long Eigen::internal::first_aligned<16, int, long>(int const*, long)
  0.00     13.49     0.00     1252     0.00     0.00  Eigen::internal::scalar_constant_op<int>::scalar_constant_op(int const&)
  0.00     13.49     0.00     1252     0.00     0.00  void Eigen::internal::unaligned_dense_assignment_loop<false>::run<Eigen::internal::generic_dense_assignment_kernel<Eigen::internal::evaluator<Eigen::Map<Eigen::Matrix<int, -1, 1, 0, -1, 1>, 0, Eigen::Stride<0, 0> > >, Eigen::internal::evaluator<Eigen::CwiseNullaryOp<Eigen::internal::scalar_constant_op<int>, Eigen::Matrix<int, -1, 1, 0, -1, 1> > >, Eigen::internal::assign_op<int, int>, 0> >(Eigen::internal::generic_dense_assignment_kernel<Eigen::internal::evaluator<Eigen::Map<Eigen::Matrix<int, -1, 1, 0, -1, 1>, 0, Eigen::Stride<0, 0> > >, Eigen::internal::evaluator<Eigen::CwiseNullaryOp<Eigen::internal::scalar_constant_op<int>, Eigen::Matrix<int, -1, 1, 0, -1, 1> > >, Eigen::internal::assign_op<int, int>, 0>&, long, long)
  0.00     13.49     0.00     1252     0.00     0.00  Eigen::internal::unpacket_traits<long long __vector(2)>::type Eigen::internal::pfirst<long long __vector(2)>(long long __vector(2) const&)
  0.00     13.49     0.00     1252     0.00     0.00  Eigen::Map<Eigen::Matrix<int, -1, 1, 0, -1, 1> const, 0, Eigen::Stride<0, 0> >::innerStride() const
  0.00     13.49     0.00     1252     0.00     0.00  Eigen::Map<Eigen::Matrix<int, -1, 1, 0, -1, 1>, 0, Eigen::Stride<0, 0> >::innerStride() const
  0.00     13.49     0.00     1252     0.00     0.00  Eigen::Stride<0, 0>::inner() const
  0.00     13.49     0.00     1252     0.00     0.00  Eigen::Stride<0, 0>::outer() const
  0.00     13.49     0.00     1252     0.00     0.00  Eigen::MapBase<Eigen::Map<Eigen::Matrix<int, -1, 1, 0, -1, 1> const, 0, Eigen::Stride<0, 0> >, 0>::data() const
  0.00     13.49     0.00     1252     0.00     0.00  Eigen::EigenBase<Eigen::Map<Eigen::Matrix<int, -1, 1, 0, -1, 1>, 0, Eigen::Stride<0, 0> > >::size() const
  0.00     13.49     0.00      636     0.00     0.00  Eigen::DenseStorage<int, -1, -1, 1, 0>::~DenseStorage()
  0.00     13.49     0.00      636     0.00     0.00  Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >::~PlainObjectBase()
  0.00     13.49     0.00      636     0.00     0.00  Eigen::Matrix<int, -1, 1, 0, -1, 1>::~Matrix()
  0.00     13.49     0.00      626     0.00     0.00  Eigen::MatrixBase<Eigen::CwiseNullaryOp<Eigen::internal::scalar_constant_op<int>, Eigen::Matrix<int, -1, 1, 0, -1, 1> > >::MatrixBase()
  0.00     13.49     0.00      626     0.00     0.00  Eigen::MatrixBase<Eigen::Map<Eigen::Matrix<int, -1, 1, 0, -1, 1> const, 0, Eigen::Stride<0, 0> > >::MatrixBase()
  0.00     13.49     0.00      626     0.00     0.00  Eigen::MatrixBase<Eigen::Map<Eigen::Matrix<int, -1, 1, 0, -1, 1>, 0, Eigen::Stride<0, 0> > >::MatrixBase()
  0.00     13.49     0.00      626     0.00     0.02  Eigen::Map<Eigen::Matrix<int, -1, 1, 0, -1, 1>, 0, Eigen::Stride<0, 0> >& Eigen::MatrixBase<Eigen::Map<Eigen::Matrix<int, -1, 1, 0, -1, 1>, 0, Eigen::Stride<0, 0> > >::operator=<Eigen::CwiseNullaryOp<Eigen::internal::scalar_constant_op<int>, Eigen::Matrix<int, -1, 1, 0, -1, 1> > >(Eigen::DenseBase<Eigen::CwiseNullaryOp<Eigen::internal::scalar_constant_op<int>, Eigen::Matrix<int, -1, 1, 0, -1, 1> > > const&)
  0.00     13.49     0.00      626     0.00     0.00  Eigen::SparseMatrix<double, 1, int>::check_template_parameters()
  0.00     13.49     0.00      626     0.00     0.00  Eigen::SparseMatrix<double, 1, int>::swap(Eigen::SparseMatrix<double, 1, int>&)
  0.00     13.49     0.00      626     0.00     0.00  Eigen::SparseMatrix<double, 1, int>::resize(long, long)
  0.00     13.49     0.00      626     0.00     0.00  Eigen::SparseMatrix<double, 1, int>::SparseMatrix(long, long)
  0.00     13.49     0.00      626     0.00     0.00  Eigen::SparseMatrix<double, 1, int>::~SparseMatrix()
  0.00     13.49     0.00      626     0.00     0.11  Eigen::SparseMatrix<double, 1, int>& Eigen::SparseMatrix<double, 1, int>::operator=<Eigen::SparseMatrix<double, 0, int> >(Eigen::SparseMatrixBase<Eigen::SparseMatrix<double, 0, int> > const&)
  0.00     13.49     0.00      626     0.00     0.00  Eigen::CwiseNullaryOp<Eigen::internal::scalar_constant_op<int>, Eigen::Array<int, -1, 1, 0, -1, 1> >::CwiseNullaryOp(long, long, Eigen::internal::scalar_constant_op<int> const&)
  0.00     13.49     0.00      626     0.00     0.00  Eigen::CwiseNullaryOp<Eigen::internal::scalar_constant_op<int>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >::CwiseNullaryOp(long, long, Eigen::internal::scalar_constant_op<int> const&)
  0.00     13.49     0.00      626     0.00     0.00  Eigen::SparseMatrixBase<Eigen::SparseMatrix<double, 1, int> >::SparseMatrixBase()
  0.00     13.49     0.00      626     0.00     0.00  Eigen::SparseCompressedBase<Eigen::SparseMatrix<double, 1, int> >::SparseCompressedBase()
  0.00     13.49     0.00      626     0.00     0.00  Eigen::Map<Eigen::Matrix<int, -1, 1, 0, -1, 1> const, 0, Eigen::Stride<0, 0> >::cast_to_pointer_type(int const*)
  0.00     13.49     0.00      626     0.00     0.00  Eigen::Map<Eigen::Matrix<int, -1, 1, 0, -1, 1> const, 0, Eigen::Stride<0, 0> >::Map(int const*, long, Eigen::Stride<0, 0> const&)
  0.00     13.49     0.00      626     0.00     0.00  Eigen::Map<Eigen::Matrix<int, -1, 1, 0, -1, 1>, 0, Eigen::Stride<0, 0> >::cast_to_pointer_type(int*)
  0.00     13.49     0.00      626     0.00     0.00  Eigen::Map<Eigen::Matrix<int, -1, 1, 0, -1, 1>, 0, Eigen::Stride<0, 0> >::Map(int*, long, Eigen::Stride<0, 0> const&)
  0.00     13.49     0.00      626     0.00     0.00  Eigen::MapBase<Eigen::Map<Eigen::Matrix<int, -1, 1, 0, -1, 1> const, 0, Eigen::Stride<0, 0> >, 0>::MapBase(int const*, long)
  0.00     13.49     0.00      626     0.00     0.00  Eigen::MapBase<Eigen::Map<Eigen::Matrix<int, -1, 1, 0, -1, 1>, 0, Eigen::Stride<0, 0> >, 0>::MapBase(int*, long)
  0.00     13.49     0.00      626     0.00     0.00  Eigen::MapBase<Eigen::Map<Eigen::Matrix<int, -1, 1, 0, -1, 1>, 0, Eigen::Stride<0, 0> >, 1>::data()
  0.00     13.49     0.00      626     0.00     0.00  Eigen::MapBase<Eigen::Map<Eigen::Matrix<int, -1, 1, 0, -1, 1>, 0, Eigen::Stride<0, 0> >, 1>::MapBase(int*, long)
  0.00     13.49     0.00      626     0.00     0.02  Eigen::internal::Assignment<Eigen::Map<Eigen::Matrix<int, -1, 1, 0, -1, 1>, 0, Eigen::Stride<0, 0> >, Eigen::CwiseNullaryOp<Eigen::internal::scalar_constant_op<int>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >, Eigen::internal::assign_op<int, int>, Eigen::internal::Dense2Dense, void>::run(Eigen::Map<Eigen::Matrix<int, -1, 1, 0, -1, 1>, 0, Eigen::Stride<0, 0> >&, Eigen::CwiseNullaryOp<Eigen::internal::scalar_constant_op<int>, Eigen::Matrix<int, -1, 1, 0, -1, 1> > const&, Eigen::internal::assign_op<int, int> const&)
  0.00     13.49     0.00      626     0.00     0.00  Eigen::internal::redux_impl<Eigen::internal::scalar_sum_op<int, int>, Eigen::internal::redux_evaluator<Eigen::Map<Eigen::Matrix<int, -1, 1, 0, -1, 1> const, 0, Eigen::Stride<0, 0> > >, 3, 0>::run(Eigen::internal::redux_evaluator<Eigen::Map<Eigen::Matrix<int, -1, 1, 0, -1, 1> const, 0, Eigen::Stride<0, 0> > > const&, Eigen::internal::scalar_sum_op<int, int> const&)
  0.00     13.49     0.00      626     0.00     0.00  Eigen::internal::scalar_sum_op<int, int>::scalar_sum_op()
  0.00     13.49     0.00      626     0.00     0.00  Eigen::internal::evaluator_base<Eigen::SparseMatrix<double, 0, int> >::evaluator_base()
  0.00     13.49     0.00      626     0.00     0.00  Eigen::internal::evaluator_base<Eigen::SparseMatrix<double, 0, int> >::~evaluator_base()
  0.00     13.49     0.00      626     0.00     0.00  Eigen::internal::evaluator_base<Eigen::CwiseNullaryOp<Eigen::internal::scalar_constant_op<int>, Eigen::Matrix<int, -1, 1, 0, -1, 1> > >::evaluator_base()
  0.00     13.49     0.00      626     0.00     0.00  Eigen::internal::evaluator_base<Eigen::CwiseNullaryOp<Eigen::internal::scalar_constant_op<int>, Eigen::Matrix<int, -1, 1, 0, -1, 1> > >::~evaluator_base()
  0.00     13.49     0.00      626     0.00     0.00  Eigen::internal::evaluator_base<Eigen::Map<Eigen::Matrix<int, -1, 1, 0, -1, 1> const, 0, Eigen::Stride<0, 0> > >::evaluator_base()
  0.00     13.49     0.00      626     0.00     0.00  Eigen::internal::evaluator_base<Eigen::Map<Eigen::Matrix<int, -1, 1, 0, -1, 1> const, 0, Eigen::Stride<0, 0> > >::~evaluator_base()
  0.00     13.49     0.00      626     0.00     0.00  Eigen::internal::evaluator_base<Eigen::Map<Eigen::Matrix<int, -1, 1, 0, -1, 1>, 0, Eigen::Stride<0, 0> > >::evaluator_base()
  0.00     13.49     0.00      626     0.00     0.00  Eigen::internal::evaluator_base<Eigen::Map<Eigen::Matrix<int, -1, 1, 0, -1, 1>, 0, Eigen::Stride<0, 0> > >::~evaluator_base()
  0.00     13.49     0.00      626     0.00     0.02  void Eigen::internal::call_assignment<Eigen::Map<Eigen::Matrix<int, -1, 1, 0, -1, 1>, 0, Eigen::Stride<0, 0> >, Eigen::CwiseNullaryOp<Eigen::internal::scalar_constant_op<int>, Eigen::Matrix<int, -1, 1, 0, -1, 1> > >(Eigen::Map<Eigen::Matrix<int, -1, 1, 0, -1, 1>, 0, Eigen::Stride<0, 0> >&, Eigen::CwiseNullaryOp<Eigen::internal::scalar_constant_op<int>, Eigen::Matrix<int, -1, 1, 0, -1, 1> > const&)
  0.00     13.49     0.00      626     0.00     0.02  void Eigen::internal::call_assignment<Eigen::Map<Eigen::Matrix<int, -1, 1, 0, -1, 1>, 0, Eigen::Stride<0, 0> >, Eigen::CwiseNullaryOp<Eigen::internal::scalar_constant_op<int>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >, Eigen::internal::assign_op<int, int> >(Eigen::Map<Eigen::Matrix<int, -1, 1, 0, -1, 1>, 0, Eigen::Stride<0, 0> >&, Eigen::CwiseNullaryOp<Eigen::internal::scalar_constant_op<int>, Eigen::Matrix<int, -1, 1, 0, -1, 1> > const&, Eigen::internal::assign_op<int, int> const&, Eigen::internal::enable_if<!Eigen::internal::evaluator_assume_aliasing<Eigen::CwiseNullaryOp<Eigen::internal::scalar_constant_op<int>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >, Eigen::internal::evaluator_traits<Eigen::CwiseNullaryOp<Eigen::internal::scalar_constant_op<int>, Eigen::Matrix<int, -1, 1, 0, -1, 1> > >::Shape>::value, void*>::type)
  0.00     13.49     0.00      626     0.00     0.00  Eigen::internal::redux_evaluator<Eigen::Map<Eigen::Matrix<int, -1, 1, 0, -1, 1> const, 0, Eigen::Stride<0, 0> > >::redux_evaluator(Eigen::Map<Eigen::Matrix<int, -1, 1, 0, -1, 1> const, 0, Eigen::Stride<0, 0> > const&)
  0.00     13.49     0.00      626     0.00     0.00  Eigen::internal::redux_evaluator<Eigen::Map<Eigen::Matrix<int, -1, 1, 0, -1, 1> const, 0, Eigen::Stride<0, 0> > >::~redux_evaluator()
  0.00     13.49     0.00      626     0.00     0.00  Eigen::internal::CompressedStorage<double, int>::swap(Eigen::internal::CompressedStorage<double, int>&)
  0.00     13.49     0.00      626     0.00     0.00  Eigen::internal::mapbase_evaluator<Eigen::Map<Eigen::Matrix<int, -1, 1, 0, -1, 1> const, 0, Eigen::Stride<0, 0> >, Eigen::Matrix<int, -1, 1, 0, -1, 1> const>::mapbase_evaluator(Eigen::Map<Eigen::Matrix<int, -1, 1, 0, -1, 1> const, 0, Eigen::Stride<0, 0> > const&)
  0.00     13.49     0.00      626     0.00     0.00  Eigen::internal::mapbase_evaluator<Eigen::Map<Eigen::Matrix<int, -1, 1, 0, -1, 1> const, 0, Eigen::Stride<0, 0> >, Eigen::Matrix<int, -1, 1, 0, -1, 1> const>::~mapbase_evaluator()
  0.00     13.49     0.00      626     0.00     0.00  Eigen::internal::mapbase_evaluator<Eigen::Map<Eigen::Matrix<int, -1, 1, 0, -1, 1>, 0, Eigen::Stride<0, 0> >, Eigen::Matrix<int, -1, 1, 0, -1, 1> >::mapbase_evaluator(Eigen::Map<Eigen::Matrix<int, -1, 1, 0, -1, 1>, 0, Eigen::Stride<0, 0> > const&)
  0.00     13.49     0.00      626     0.00     0.00  Eigen::internal::mapbase_evaluator<Eigen::Map<Eigen::Matrix<int, -1, 1, 0, -1, 1>, 0, Eigen::Stride<0, 0> >, Eigen::Matrix<int, -1, 1, 0, -1, 1> >::~mapbase_evaluator()
  0.00     13.49     0.00      626     0.00     0.00  void Eigen::internal::resize_if_allowed<Eigen::Map<Eigen::Matrix<int, -1, 1, 0, -1, 1>, 0, Eigen::Stride<0, 0> >, Eigen::CwiseNullaryOp<Eigen::internal::scalar_constant_op<int>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >, int, int>(Eigen::Map<Eigen::Matrix<int, -1, 1, 0, -1, 1>, 0, Eigen::Stride<0, 0> >&, Eigen::CwiseNullaryOp<Eigen::internal::scalar_constant_op<int>, Eigen::Matrix<int, -1, 1, 0, -1, 1> > const&, Eigen::internal::assign_op<int, int> const&)
  0.00     13.49     0.00      626     0.00     0.00  void Eigen::internal::check_for_aliasing<Eigen::Map<Eigen::Matrix<int, -1, 1, 0, -1, 1>, 0, Eigen::Stride<0, 0> >, Eigen::CwiseNullaryOp<Eigen::internal::scalar_constant_op<int>, Eigen::Matrix<int, -1, 1, 0, -1, 1> > >(Eigen::Map<Eigen::Matrix<int, -1, 1, 0, -1, 1>, 0, Eigen::Stride<0, 0> > const&, Eigen::CwiseNullaryOp<Eigen::internal::scalar_constant_op<int>, Eigen::Matrix<int, -1, 1, 0, -1, 1> > const&)
  0.00     13.49     0.00      626     0.00     0.00  Eigen::internal::first_aligned_impl<16, Eigen::Map<Eigen::Matrix<int, -1, 1, 0, -1, 1> const, 0, Eigen::Stride<0, 0> >, false>::run(Eigen::Map<Eigen::Matrix<int, -1, 1, 0, -1, 1> const, 0, Eigen::Stride<0, 0> > const&)
  0.00     13.49     0.00      626     0.00     0.02  Eigen::internal::dense_assignment_loop<Eigen::internal::generic_dense_assignment_kernel<Eigen::internal::evaluator<Eigen::Map<Eigen::Matrix<int, -1, 1, 0, -1, 1>, 0, Eigen::Stride<0, 0> > >, Eigen::internal::evaluator<Eigen::CwiseNullaryOp<Eigen::internal::scalar_constant_op<int>, Eigen::Matrix<int, -1, 1, 0, -1, 1> > >, Eigen::internal::assign_op<int, int>, 0>, 3, 0>::run(Eigen::internal::generic_dense_assignment_kernel<Eigen::internal::evaluator<Eigen::Map<Eigen::Matrix<int, -1, 1, 0, -1, 1>, 0, Eigen::Stride<0, 0> > >, Eigen::internal::evaluator<Eigen::CwiseNullaryOp<Eigen::internal::scalar_constant_op<int>, Eigen::Matrix<int, -1, 1, 0, -1, 1> > >, Eigen::internal::assign_op<int, int>, 0>&)
  0.00     13.49     0.00      626     0.00     0.02  void Eigen::internal::call_assignment_no_alias<Eigen::Map<Eigen::Matrix<int, -1, 1, 0, -1, 1>, 0, Eigen::Stride<0, 0> >, Eigen::CwiseNullaryOp<Eigen::internal::scalar_constant_op<int>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >, Eigen::internal::assign_op<int, int> >(Eigen::Map<Eigen::Matrix<int, -1, 1, 0, -1, 1>, 0, Eigen::Stride<0, 0> >&, Eigen::CwiseNullaryOp<Eigen::internal::scalar_constant_op<int>, Eigen::Matrix<int, -1, 1, 0, -1, 1> > const&, Eigen::internal::assign_op<int, int> const&)
  0.00     13.49     0.00      626     0.00     0.02  void Eigen::internal::call_dense_assignment_loop<Eigen::Map<Eigen::Matrix<int, -1, 1, 0, -1, 1>, 0, Eigen::Stride<0, 0> >, Eigen::CwiseNullaryOp<Eigen::internal::scalar_constant_op<int>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >, Eigen::internal::assign_op<int, int> >(Eigen::Map<Eigen::Matrix<int, -1, 1, 0, -1, 1>, 0, Eigen::Stride<0, 0> >&, Eigen::CwiseNullaryOp<Eigen::internal::scalar_constant_op<int>, Eigen::Matrix<int, -1, 1, 0, -1, 1> > const&, Eigen::internal::assign_op<int, int> const&)
  0.00     13.49     0.00      626     0.00     0.00  Eigen::internal::checkTransposeAliasing_impl<Eigen::Map<Eigen::Matrix<int, -1, 1, 0, -1, 1>, 0, Eigen::Stride<0, 0> >, Eigen::CwiseNullaryOp<Eigen::internal::scalar_constant_op<int>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >, false>::run(Eigen::Map<Eigen::Matrix<int, -1, 1, 0, -1, 1>, 0, Eigen::Stride<0, 0> > const&, Eigen::CwiseNullaryOp<Eigen::internal::scalar_constant_op<int>, Eigen::Matrix<int, -1, 1, 0, -1, 1> > const&)
  0.00     13.49     0.00      626     0.00     0.00  Eigen::internal::generic_dense_assignment_kernel<Eigen::internal::evaluator<Eigen::Map<Eigen::Matrix<int, -1, 1, 0, -1, 1>, 0, Eigen::Stride<0, 0> > >, Eigen::internal::evaluator<Eigen::CwiseNullaryOp<Eigen::internal::scalar_constant_op<int>, Eigen::Matrix<int, -1, 1, 0, -1, 1> > >, Eigen::internal::assign_op<int, int>, 0>::assignCoeff(long)
  0.00     13.49     0.00      626     0.00     0.00  Eigen::internal::generic_dense_assignment_kernel<Eigen::internal::evaluator<Eigen::Map<Eigen::Matrix<int, -1, 1, 0, -1, 1>, 0, Eigen::Stride<0, 0> > >, Eigen::internal::evaluator<Eigen::CwiseNullaryOp<Eigen::internal::scalar_constant_op<int>, Eigen::Matrix<int, -1, 1, 0, -1, 1> > >, Eigen::internal::assign_op<int, int>, 0>::generic_dense_assignment_kernel(Eigen::internal::evaluator<Eigen::Map<Eigen::Matrix<int, -1, 1, 0, -1, 1>, 0, Eigen::Stride<0, 0> > >&, Eigen::internal::evaluator<Eigen::CwiseNullaryOp<Eigen::internal::scalar_constant_op<int>, Eigen::Matrix<int, -1, 1, 0, -1, 1> > > const&, Eigen::internal::assign_op<int, int> const&, Eigen::Map<Eigen::Matrix<int, -1, 1, 0, -1, 1>, 0, Eigen::Stride<0, 0> >&)
  0.00     13.49     0.00      626     0.00     0.00  Eigen::internal::unpacket_traits<long long __vector(2)>::type Eigen::internal::predux<long long __vector(2)>(long long __vector(2) const&)
  0.00     13.49     0.00      626     0.00     0.00  Eigen::internal::assign_op<int, int>::assign_op()
  0.00     13.49     0.00      626     0.00     0.00  Eigen::internal::evaluator<Eigen::SparseMatrix<double, 0, int> >::evaluator(Eigen::SparseMatrix<double, 0, int> const&)
  0.00     13.49     0.00      626     0.00     0.00  Eigen::internal::evaluator<Eigen::SparseMatrix<double, 0, int> >::~evaluator()
  0.00     13.49     0.00      626     0.00     0.00  Eigen::internal::evaluator<Eigen::CwiseNullaryOp<Eigen::internal::scalar_constant_op<int>, Eigen::Matrix<int, -1, 1, 0, -1, 1> > >::evaluator(Eigen::CwiseNullaryOp<Eigen::internal::scalar_constant_op<int>, Eigen::Matrix<int, -1, 1, 0, -1, 1> > const&)
  0.00     13.49     0.00      626     0.00     0.00  Eigen::internal::evaluator<Eigen::CwiseNullaryOp<Eigen::internal::scalar_constant_op<int>, Eigen::Matrix<int, -1, 1, 0, -1, 1> > >::~evaluator()
  0.00     13.49     0.00      626     0.00     0.00  Eigen::internal::evaluator<Eigen::SparseCompressedBase<Eigen::SparseMatrix<double, 0, int> > >::evaluator(Eigen::SparseMatrix<double, 0, int> const&)
  0.00     13.49     0.00      626     0.00     0.00  Eigen::internal::evaluator<Eigen::SparseCompressedBase<Eigen::SparseMatrix<double, 0, int> > >::~evaluator()
  0.00     13.49     0.00      626     0.00     0.00  Eigen::internal::evaluator<Eigen::Map<Eigen::Matrix<int, -1, 1, 0, -1, 1> const, 0, Eigen::Stride<0, 0> > >::evaluator(Eigen::Map<Eigen::Matrix<int, -1, 1, 0, -1, 1> const, 0, Eigen::Stride<0, 0> > const&)
  0.00     13.49     0.00      626     0.00     0.00  Eigen::internal::evaluator<Eigen::Map<Eigen::Matrix<int, -1, 1, 0, -1, 1> const, 0, Eigen::Stride<0, 0> > >::~evaluator()
  0.00     13.49     0.00      626     0.00     0.00  Eigen::internal::evaluator<Eigen::Map<Eigen::Matrix<int, -1, 1, 0, -1, 1>, 0, Eigen::Stride<0, 0> > >::evaluator(Eigen::Map<Eigen::Matrix<int, -1, 1, 0, -1, 1>, 0, Eigen::Stride<0, 0> > const&)
  0.00     13.49     0.00      626     0.00     0.00  Eigen::internal::evaluator<Eigen::Map<Eigen::Matrix<int, -1, 1, 0, -1, 1>, 0, Eigen::Stride<0, 0> > >::~evaluator()
  0.00     13.49     0.00      626     0.00     0.00  long Eigen::internal::first_aligned<16, Eigen::Map<Eigen::Matrix<int, -1, 1, 0, -1, 1> const, 0, Eigen::Stride<0, 0> > >(Eigen::DenseBase<Eigen::Map<Eigen::Matrix<int, -1, 1, 0, -1, 1> const, 0, Eigen::Stride<0, 0> > > const&)
  0.00     13.49     0.00      626     0.00     0.00  long Eigen::internal::first_default_aligned<Eigen::Map<Eigen::Matrix<int, -1, 1, 0, -1, 1> const, 0, Eigen::Stride<0, 0> > >(Eigen::DenseBase<Eigen::Map<Eigen::Matrix<int, -1, 1, 0, -1, 1> const, 0, Eigen::Stride<0, 0> > > const&)
  0.00     13.49     0.00      626     0.00     0.00  Eigen::ArrayBase<Eigen::CwiseNullaryOp<Eigen::internal::scalar_constant_op<int>, Eigen::Array<int, -1, 1, 0, -1, 1> > >::ArrayBase()
  0.00     13.49     0.00      626     0.00     0.00  Eigen::DenseBase<Eigen::CwiseNullaryOp<Eigen::internal::scalar_constant_op<int>, Eigen::Array<int, -1, 1, 0, -1, 1> > >::DenseBase()
  0.00     13.49     0.00      626     0.00     0.00  Eigen::DenseBase<Eigen::CwiseNullaryOp<Eigen::internal::scalar_constant_op<int>, Eigen::Matrix<int, -1, 1, 0, -1, 1> > >::DenseBase()
  0.00     13.49     0.00      626     0.00     0.00  Eigen::DenseBase<Eigen::Map<Eigen::Matrix<int, -1, 1, 0, -1, 1> const, 0, Eigen::Stride<0, 0> > >::DenseBase()
  0.00     13.49     0.00      626     0.00     0.00  Eigen::CwiseNullaryOp<Eigen::internal::scalar_constant_op<int>, Eigen::Matrix<int, -1, 1, 0, -1, 1> > const Eigen::DenseBase<Eigen::Map<Eigen::Matrix<int, -1, 1, 0, -1, 1>, 0, Eigen::Stride<0, 0> > >::NullaryExpr<Eigen::internal::scalar_constant_op<int> >(long, long, Eigen::internal::scalar_constant_op<int> const&)
  0.00     13.49     0.00      626     0.00     0.02  Eigen::DenseBase<Eigen::Map<Eigen::Matrix<int, -1, 1, 0, -1, 1>, 0, Eigen::Stride<0, 0> > >::setConstant(int const&)
  0.00     13.49     0.00      626     0.00     0.02  Eigen::DenseBase<Eigen::Map<Eigen::Matrix<int, -1, 1, 0, -1, 1>, 0, Eigen::Stride<0, 0> > >::setZero()
  0.00     13.49     0.00      626     0.00     0.00  Eigen::DenseBase<Eigen::Map<Eigen::Matrix<int, -1, 1, 0, -1, 1>, 0, Eigen::Stride<0, 0> > >::Constant(long, long, int const&)
  0.00     13.49     0.00      626     0.00     0.00  Eigen::DenseBase<Eigen::Map<Eigen::Matrix<int, -1, 1, 0, -1, 1>, 0, Eigen::Stride<0, 0> > >::DenseBase()
  0.00     13.49     0.00      626     0.00     0.00  Eigen::CwiseNullaryOp<Eigen::internal::scalar_constant_op<int>, Eigen::Array<int, -1, 1, 0, -1, 1> > const Eigen::DenseBase<Eigen::Array<int, -1, 1, 0, -1, 1> >::NullaryExpr<Eigen::internal::scalar_constant_op<int> >(long, Eigen::internal::scalar_constant_op<int> const&)
  0.00     13.49     0.00      626     0.00     0.00  Eigen::DenseBase<Eigen::Array<int, -1, 1, 0, -1, 1> >::Constant(long, int const&)
  0.00     13.49     0.00      626     0.00     0.00  Eigen::CwiseNullaryOp<Eigen::internal::scalar_constant_op<int>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >::cols() const
  0.00     13.49     0.00      626     0.00     0.00  Eigen::CwiseNullaryOp<Eigen::internal::scalar_constant_op<int>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >::rows() const
  0.00     13.49     0.00      626     0.00     0.00  Eigen::CwiseNullaryOp<Eigen::internal::scalar_constant_op<int>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >::functor() const
  0.00     13.49     0.00      626     0.00     0.00  Eigen::SparseMatrixBase<Eigen::SparseMatrix<double, 0, int> >::cols() const
  0.00     13.49     0.00      626     0.00     0.00  Eigen::SparseMatrixBase<Eigen::SparseMatrix<double, 0, int> >::rows() const
  0.00     13.49     0.00      626     0.00     0.00  Eigen::SparseCompressedBase<Eigen::SparseMatrix<double, 0, int> >::innerNonZeros() const
  0.00     13.49     0.00      626     0.00     0.00  Eigen::Map<Eigen::Matrix<int, -1, 1, 0, -1, 1> const, 0, Eigen::Stride<0, 0> >::outerStride() const
  0.00     13.49     0.00      626     0.00     0.00  Eigen::Map<Eigen::Matrix<int, -1, 1, 0, -1, 1>, 0, Eigen::Stride<0, 0> >::outerStride() const
  0.00     13.49     0.00      626     0.00     0.00  void Eigen::MapBase<Eigen::Map<Eigen::Matrix<int, -1, 1, 0, -1, 1> const, 0, Eigen::Stride<0, 0> >, 0>::checkSanity<Eigen::Map<Eigen::Matrix<int, -1, 1, 0, -1, 1> const, 0, Eigen::Stride<0, 0> > >(Eigen::internal::enable_if<Eigen::internal::traits<Eigen::Map<Eigen::Matrix<int, -1, 1, 0, -1, 1> const, 0, Eigen::Stride<0, 0> > >::Alignment==(0), void*>::type) const
  0.00     13.49     0.00      626     0.00     0.00  void Eigen::MapBase<Eigen::Map<Eigen::Matrix<int, -1, 1, 0, -1, 1>, 0, Eigen::Stride<0, 0> >, 0>::checkSanity<Eigen::Map<Eigen::Matrix<int, -1, 1, 0, -1, 1>, 0, Eigen::Stride<0, 0> > >(Eigen::internal::enable_if<Eigen::internal::traits<Eigen::Map<Eigen::Matrix<int, -1, 1, 0, -1, 1>, 0, Eigen::Stride<0, 0> > >::Alignment==(0), void*>::type) const
  0.00     13.49     0.00      626     0.00     0.00  Eigen::MapBase<Eigen::Map<Eigen::Matrix<int, -1, 1, 0, -1, 1>, 0, Eigen::Stride<0, 0> >, 1>::data() const
  0.00     13.49     0.00      626     0.00     0.00  int const Eigen::internal::scalar_sum_op<int, int>::predux<long long __vector(2)>(long long __vector(2) const&) const
  0.00     13.49     0.00      626     0.00     0.00  Eigen::internal::redux_evaluator<Eigen::Map<Eigen::Matrix<int, -1, 1, 0, -1, 1> const, 0, Eigen::Stride<0, 0> > >::nestedExpression() const
  0.00     13.49     0.00      626     0.00     0.00  Eigen::internal::redux_evaluator<Eigen::Map<Eigen::Matrix<int, -1, 1, 0, -1, 1> const, 0, Eigen::Stride<0, 0> > >::size() const
  0.00     13.49     0.00      626     0.00     0.00  Eigen::internal::generic_dense_assignment_kernel<Eigen::internal::evaluator<Eigen::Map<Eigen::Matrix<int, -1, 1, 0, -1, 1>, 0, Eigen::Stride<0, 0> > >, Eigen::internal::evaluator<Eigen::CwiseNullaryOp<Eigen::internal::scalar_constant_op<int>, Eigen::Matrix<int, -1, 1, 0, -1, 1> > >, Eigen::internal::assign_op<int, int>, 0>::dstDataPtr() const
  0.00     13.49     0.00      626     0.00     0.00  Eigen::internal::generic_dense_assignment_kernel<Eigen::internal::evaluator<Eigen::Map<Eigen::Matrix<int, -1, 1, 0, -1, 1>, 0, Eigen::Stride<0, 0> > >, Eigen::internal::evaluator<Eigen::CwiseNullaryOp<Eigen::internal::scalar_constant_op<int>, Eigen::Matrix<int, -1, 1, 0, -1, 1> > >, Eigen::internal::assign_op<int, int>, 0>::size() const
  0.00     13.49     0.00      626     0.00     0.00  Eigen::internal::assign_op<int, int>::assignCoeff(int&, int const&) const
  0.00     13.49     0.00      626     0.00     0.00  int Eigen::internal::evaluator<Eigen::CwiseNullaryOp<Eigen::internal::scalar_constant_op<int>, Eigen::Matrix<int, -1, 1, 0, -1, 1> > >::coeff<long>(long) const
  0.00     13.49     0.00      626     0.00     0.00  Eigen::DenseBase<Eigen::Map<Eigen::Matrix<int, -1, 1, 0, -1, 1> const, 0, Eigen::Stride<0, 0> > >::sum() const
  0.00     13.49     0.00      626     0.00     0.00  int Eigen::DenseBase<Eigen::Map<Eigen::Matrix<int, -1, 1, 0, -1, 1> const, 0, Eigen::Stride<0, 0> > >::redux<Eigen::internal::scalar_sum_op<int, int> >(Eigen::internal::scalar_sum_op<int, int> const&) const
  0.00     13.49     0.00      626     0.00     0.00  Eigen::EigenBase<Eigen::CwiseNullaryOp<Eigen::internal::scalar_constant_op<int>, Eigen::Matrix<int, -1, 1, 0, -1, 1> > >::derived() const
  0.00     13.49     0.00      626     0.00     0.00  Eigen::EigenBase<Eigen::Map<Eigen::Matrix<int, -1, 1, 0, -1, 1>, 0, Eigen::Stride<0, 0> > >::const_cast_derived() const
  0.00     13.49     0.00       10     0.00     0.00  Eigen::Matrix<int, -1, 1, 0, -1, 1>::Matrix()
  0.00     13.49     0.00        9     0.00     0.00  __gnu_cxx::new_allocator<double>::~new_allocator()
  0.00     13.49     0.00        9     0.00     0.00  std::allocator<double>::~allocator()
  0.00     13.49     0.00        9     0.00     0.00  void std::_Destroy_aux<true>::__destroy<double*>(double*, double*)
  0.00     13.49     0.00        9     0.00     0.00  std::_Vector_base<double, std::allocator<double> >::_Vector_impl::~_Vector_impl()
  0.00     13.49     0.00        9     0.00     0.00  std::_Vector_base<double, std::allocator<double> >::_M_deallocate(double*, unsigned long)
  0.00     13.49     0.00        9     0.00     0.00  std::_Vector_base<double, std::allocator<double> >::_M_get_Tp_allocator()
  0.00     13.49     0.00        9     0.00     0.00  std::_Vector_base<double, std::allocator<double> >::~_Vector_base()
  0.00     13.49     0.00        9     0.00     0.00  std::vector<double, std::allocator<double> >::~vector()
  0.00     13.49     0.00        9     0.00     0.00  void std::_Destroy<double*>(double*, double*)
  0.00     13.49     0.00        9     0.00     0.00  void std::_Destroy<double*, double>(double*, double*, std::allocator<double>&)
  0.00     13.49     0.00        6     0.00     0.00  __gnu_cxx::new_allocator<Cantera::FlowDevice*>::~new_allocator()
  0.00     13.49     0.00        6     0.00     0.00  __gnu_cxx::new_allocator<double>::deallocate(double*, unsigned long)
  0.00     13.49     0.00        6     0.00     0.00  std::allocator<Cantera::FlowDevice*>::~allocator()
  0.00     13.49     0.00        6     0.00     0.00  void std::_Destroy_aux<true>::__destroy<Cantera::FlowDevice**>(Cantera::FlowDevice**, Cantera::FlowDevice**)
  0.00     13.49     0.00        6     0.00     0.00  std::_Vector_base<Cantera::FlowDevice*, std::allocator<Cantera::FlowDevice*> >::_Vector_impl::~_Vector_impl()
  0.00     13.49     0.00        6     0.00     0.00  std::_Vector_base<Cantera::FlowDevice*, std::allocator<Cantera::FlowDevice*> >::_M_deallocate(Cantera::FlowDevice**, unsigned long)
  0.00     13.49     0.00        6     0.00     0.00  std::_Vector_base<Cantera::FlowDevice*, std::allocator<Cantera::FlowDevice*> >::_M_get_Tp_allocator()
  0.00     13.49     0.00        6     0.00     0.00  std::_Vector_base<Cantera::FlowDevice*, std::allocator<Cantera::FlowDevice*> >::~_Vector_base()
  0.00     13.49     0.00        6     0.00     0.00  std::allocator_traits<std::allocator<double> >::deallocate(std::allocator<double>&, double*, unsigned long)
  0.00     13.49     0.00        6     0.00     0.00  std::vector<Cantera::FlowDevice*, std::allocator<Cantera::FlowDevice*> >::~vector()
  0.00     13.49     0.00        6     0.00     0.00  void std::_Destroy<Cantera::FlowDevice**>(Cantera::FlowDevice**, Cantera::FlowDevice**)
  0.00     13.49     0.00        6     0.00     0.00  void std::_Destroy<Cantera::FlowDevice**, Cantera::FlowDevice*>(Cantera::FlowDevice**, Cantera::FlowDevice**, std::allocator<Cantera::FlowDevice*>&)
  0.00     13.49     0.00        5     0.00     0.00  __gnu_cxx::new_allocator<unsigned long>::deallocate(unsigned long*, unsigned long)
  0.00     13.49     0.00        5     0.00     0.00  __gnu_cxx::new_allocator<unsigned long>::~new_allocator()
  0.00     13.49     0.00        5     0.00     0.00  std::allocator<unsigned long>::~allocator()
  0.00     13.49     0.00        5     0.00     0.00  void std::_Destroy_aux<true>::__destroy<unsigned long*>(unsigned long*, unsigned long*)
  0.00     13.49     0.00        5     0.00     0.00  std::_Vector_base<unsigned long, std::allocator<unsigned long> >::_Vector_impl::~_Vector_impl()
  0.00     13.49     0.00        5     0.00     0.00  std::_Vector_base<unsigned long, std::allocator<unsigned long> >::_M_deallocate(unsigned long*, unsigned long)
  0.00     13.49     0.00        5     0.00     0.00  std::_Vector_base<unsigned long, std::allocator<unsigned long> >::_M_get_Tp_allocator()
  0.00     13.49     0.00        5     0.00     0.00  std::_Vector_base<unsigned long, std::allocator<unsigned long> >::~_Vector_base()
  0.00     13.49     0.00        5     0.00     0.00  std::allocator_traits<std::allocator<unsigned long> >::deallocate(std::allocator<unsigned long>&, unsigned long*, unsigned long)
  0.00     13.49     0.00        5     0.00     0.00  std::vector<unsigned long, std::allocator<unsigned long> >::~vector()
  0.00     13.49     0.00        5     0.00     0.00  void std::_Destroy<unsigned long*>(unsigned long*, unsigned long*)
  0.00     13.49     0.00        5     0.00     0.00  void std::_Destroy<unsigned long*, unsigned long>(unsigned long*, unsigned long*, std::allocator<unsigned long>&)
  0.00     13.49     0.00        4     0.00     0.00  __gnu_cxx::new_allocator<Cantera::FlowDevice*>::deallocate(Cantera::FlowDevice**, unsigned long)
  0.00     13.49     0.00        4     0.00     0.00  std::__shared_ptr<Cantera::ThermoPhase, (__gnu_cxx::_Lock_policy)2>::get() const
  0.00     13.49     0.00        4     0.00     0.00  std::__shared_ptr_access<Cantera::ThermoPhase, (__gnu_cxx::_Lock_policy)2, false, false>::_M_get() const
  0.00     13.49     0.00        4     0.00     0.00  std::allocator_traits<std::allocator<Cantera::FlowDevice*> >::deallocate(std::allocator<Cantera::FlowDevice*>&, Cantera::FlowDevice**, unsigned long)
  0.00     13.49     0.00        3     0.00     0.00  Eigen::SparseMatrix<double, 0, int>::SparseMatrix()
  0.00     13.49     0.00        3     0.00     0.00  Cantera::ReactorBase::~ReactorBase()
  0.00     13.49     0.00        3     0.00     0.00  __gnu_cxx::new_allocator<Cantera::ReactorSurface*>::~new_allocator()
  0.00     13.49     0.00        3     0.00     0.00  __gnu_cxx::new_allocator<Cantera::WallBase*>::~new_allocator()
  0.00     13.49     0.00        3     0.00     0.00  __gnu_cxx::new_allocator<int>::~new_allocator()
  0.00     13.49     0.00        3     0.00     0.00  std::__shared_ptr_access<Cantera::ThermoPhase, (__gnu_cxx::_Lock_policy)2, false, false>::operator*() const
  0.00     13.49     0.00        3     0.00     0.00  std::allocator<Cantera::ReactorSurface*>::~allocator()
  0.00     13.49     0.00        3     0.00     0.00  std::allocator<Cantera::WallBase*>::~allocator()
  0.00     13.49     0.00        3     0.00     0.00  std::allocator<int>::~allocator()
  0.00     13.49     0.00        3     0.00     0.00  void std::_Destroy_aux<true>::__destroy<Cantera::ReactorSurface**>(Cantera::ReactorSurface**, Cantera::ReactorSurface**)
  0.00     13.49     0.00        3     0.00     0.00  void std::_Destroy_aux<true>::__destroy<Cantera::WallBase**>(Cantera::WallBase**, Cantera::WallBase**)
  0.00     13.49     0.00        3     0.00     0.00  void std::_Destroy_aux<true>::__destroy<int*>(int*, int*)
  0.00     13.49     0.00        3     0.00     0.00  std::_Vector_base<Cantera::ReactorSurface*, std::allocator<Cantera::ReactorSurface*> >::_Vector_impl::~_Vector_impl()
  0.00     13.49     0.00        3     0.00     0.00  std::_Vector_base<Cantera::ReactorSurface*, std::allocator<Cantera::ReactorSurface*> >::_M_deallocate(Cantera::ReactorSurface**, unsigned long)
  0.00     13.49     0.00        3     0.00     0.00  std::_Vector_base<Cantera::ReactorSurface*, std::allocator<Cantera::ReactorSurface*> >::_M_get_Tp_allocator()
  0.00     13.49     0.00        3     0.00     0.00  std::_Vector_base<Cantera::ReactorSurface*, std::allocator<Cantera::ReactorSurface*> >::~_Vector_base()
  0.00     13.49     0.00        3     0.00     0.00  std::_Vector_base<Cantera::WallBase*, std::allocator<Cantera::WallBase*> >::_Vector_impl::~_Vector_impl()
  0.00     13.49     0.00        3     0.00     0.00  std::_Vector_base<Cantera::WallBase*, std::allocator<Cantera::WallBase*> >::_M_deallocate(Cantera::WallBase**, unsigned long)
  0.00     13.49     0.00        3     0.00     0.00  std::_Vector_base<Cantera::WallBase*, std::allocator<Cantera::WallBase*> >::_M_get_Tp_allocator()
  0.00     13.49     0.00        3     0.00     0.00  std::_Vector_base<Cantera::WallBase*, std::allocator<Cantera::WallBase*> >::~_Vector_base()
  0.00     13.49     0.00        3     0.00     0.00  std::_Vector_base<int, std::allocator<int> >::_Vector_impl::~_Vector_impl()
  0.00     13.49     0.00        3     0.00     0.00  std::_Vector_base<int, std::allocator<int> >::_M_deallocate(int*, unsigned long)
  0.00     13.49     0.00        3     0.00     0.00  std::_Vector_base<int, std::allocator<int> >::_M_get_Tp_allocator()
  0.00     13.49     0.00        3     0.00     0.00  std::_Vector_base<int, std::allocator<int> >::~_Vector_base()
  0.00     13.49     0.00        3     0.00     0.00  std::__shared_count<(__gnu_cxx::_Lock_policy)2>::~__shared_count()
  0.00     13.49     0.00        3     0.00     0.00  std::vector<Cantera::ReactorSurface*, std::allocator<Cantera::ReactorSurface*> >::~vector()
  0.00     13.49     0.00        3     0.00     0.00  std::vector<Cantera::WallBase*, std::allocator<Cantera::WallBase*> >::~vector()
  0.00     13.49     0.00        3     0.00     0.00  std::vector<int, std::allocator<int> >::~vector()
  0.00     13.49     0.00        3     0.00     0.00  void std::_Destroy<Cantera::ReactorSurface**>(Cantera::ReactorSurface**, Cantera::ReactorSurface**)
  0.00     13.49     0.00        3     0.00     0.00  void std::_Destroy<Cantera::ReactorSurface**, Cantera::ReactorSurface*>(Cantera::ReactorSurface**, Cantera::ReactorSurface**, std::allocator<Cantera::ReactorSurface*>&)
  0.00     13.49     0.00        3     0.00     0.00  void std::_Destroy<Cantera::WallBase**>(Cantera::WallBase**, Cantera::WallBase**)
  0.00     13.49     0.00        3     0.00     0.00  void std::_Destroy<Cantera::WallBase**, Cantera::WallBase*>(Cantera::WallBase**, Cantera::WallBase**, std::allocator<Cantera::WallBase*>&)
  0.00     13.49     0.00        3     0.00     0.00  void std::_Destroy<int*>(int*, int*)
  0.00     13.49     0.00        3     0.00     0.00  void std::_Destroy<int*, int>(int*, int*, std::allocator<int>&)
  0.00     13.49     0.00        2     0.00     0.00  Eigen::MatrixBase<Eigen::Matrix<double, -1, 1, 0, -1, 1> >::MatrixBase()
  0.00     13.49     0.00        2     0.00     0.00  Eigen::DenseStorage<double, -1, -1, 1, 0>::DenseStorage()
  0.00     13.49     0.00        2     0.00     0.00  Eigen::DenseStorage<double, -1, -1, 1, 0>::~DenseStorage()
  0.00     13.49     0.00        2     0.00     0.00  Eigen::DenseStorage<int, 2, 2, 1, 0>::data()
  0.00     13.49     0.00        2     0.00     0.00  Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 1, 0, -1, 1> >::_check_template_params()
  0.00     13.49     0.00        2     0.00     0.00  Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 1, 0, -1, 1> >::PlainObjectBase()
  0.00     13.49     0.00        2     0.00     0.00  Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 1, 0, -1, 1> >::~PlainObjectBase()
  0.00     13.49     0.00        2     0.00     0.00  Eigen::PermutationMatrix<-1, -1, int>::PermutationMatrix()
  0.00     13.49     0.00        2     0.00     0.00  Eigen::PermutationMatrix<-1, -1, int>::~PermutationMatrix()
  0.00     13.49     0.00        2     0.00     0.00  Eigen::Matrix<double, -1, 1, 0, -1, 1>::Matrix()
  0.00     13.49     0.00        2     0.00     0.00  Eigen::Matrix<double, -1, 1, 0, -1, 1>::~Matrix()
  0.00     13.49     0.00        2     0.00     0.00  void Eigen::internal::conditional_aligned_delete_auto<double, true>(double*, unsigned long)
  0.00     13.49     0.00        2     0.00     0.00  Eigen::DenseBase<Eigen::Matrix<double, -1, 1, 0, -1, 1> >::DenseBase()
  0.00     13.49     0.00        2     0.00     0.00  Cantera::FlowDevice::~FlowDevice()
  0.00     13.49     0.00        2     0.00     0.00  Cantera::Reservoir::insert(Cantera::ThermoPhase&)
  0.00     13.49     0.00        2     0.00     0.00  Cantera::Reservoir::Reservoir()
  0.00     13.49     0.00        2     0.00     0.00  Cantera::Reservoir::~Reservoir()
  0.00     13.49     0.00        2     0.00     0.00  __gnu_cxx::__atomic_add(int volatile*, int)
  0.00     13.49     0.00        2     0.00     0.00  __gnu_cxx::__atomic_add_dispatch(int*, int)
  0.00     13.49     0.00        2     0.00     0.00  std::__shared_ptr<Cantera::Solution, (__gnu_cxx::_Lock_policy)2>::get() const
  0.00     13.49     0.00        2     0.00     0.00  std::__shared_ptr_access<Cantera::Solution, (__gnu_cxx::_Lock_policy)2, false, false>::_M_get() const
  0.00     13.49     0.00        2     0.00     0.00  std::__shared_ptr_access<Cantera::Solution, (__gnu_cxx::_Lock_policy)2, false, false>::operator->() const
  0.00     13.49     0.00        2     0.00     0.00  std::_Sp_counted_base<(__gnu_cxx::_Lock_policy)2>::_M_add_ref_copy()
  0.00     13.49     0.00        1     0.00     0.00  _GLOBAL__sub_I__Z20HydrogenAutoIgnitionv
  0.00     13.49     0.00        1     0.00     0.00  _GLOBAL__sub_I_main
  0.00     13.49     0.00        1     0.00     0.00  __static_initialization_and_destruction_0(int, int)
  0.00     13.49     0.00        1     0.00     0.00  __static_initialization_and_destruction_0(int, int)
  0.00     13.49     0.00        1     0.00     5.00  JetA()
  0.00     13.49     0.00        1     0.00     0.00  Eigen::DenseStorage<int, 2, 2, 1, 0>::DenseStorage()
  0.00     13.49     0.00        1     0.00     0.00  Eigen::SparseMapBase<Eigen::Map<Eigen::SparseMatrix<double, 0, int>, 0, Eigen::Stride<0, 0> >, 0>::SparseMapBase(long, long, long, int*, int*, double*, int*)
  0.00     13.49     0.00        1     0.00     0.00  Eigen::SparseMapBase<Eigen::Map<Eigen::SparseMatrix<double, 0, int>, 0, Eigen::Stride<0, 0> >, 0>::~SparseMapBase()
  0.00     13.49     0.00        1     0.00     0.00  Eigen::SparseMapBase<Eigen::Map<Eigen::SparseMatrix<double, 0, int>, 0, Eigen::Stride<0, 0> >, 1>::SparseMapBase(long, long, long, int*, int*, double*, int*)
  0.00     13.49     0.00        1     0.00     0.00  Eigen::SparseMapBase<Eigen::Map<Eigen::SparseMatrix<double, 0, int>, 0, Eigen::Stride<0, 0> >, 1>::~SparseMapBase()
  0.00     13.49     0.00        1     0.00     0.00  Eigen::PlainObjectBase<Eigen::Array<int, 2, 1, 0, 2, 1> >::_check_template_params()
  0.00     13.49     0.00        1     0.00     0.00  void Eigen::PlainObjectBase<Eigen::Array<int, 2, 1, 0, 2, 1> >::_init2<int, int>(int const&, int const&, Eigen::internal::enable_if<true, int>::type*)
  0.00     13.49     0.00        1     0.00     0.00  Eigen::PlainObjectBase<Eigen::Array<int, 2, 1, 0, 2, 1> >::PlainObjectBase()
  0.00     13.49     0.00        1     0.00     0.00  Eigen::SparseMatrixBase<Eigen::Map<Eigen::SparseMatrix<double, 0, int>, 0, Eigen::Stride<0, 0> > >::SparseMatrixBase()
  0.00     13.49     0.00        1     0.00     0.00  Eigen::SparseSolverBase<Eigen::SparseLU<Eigen::SparseMatrix<double, 0, int>, Eigen::COLAMDOrdering<int> > >::SparseSolverBase()
  0.00     13.49     0.00        1     0.00     0.00  Eigen::SparseSolverBase<Eigen::SparseLU<Eigen::SparseMatrix<double, 0, int>, Eigen::COLAMDOrdering<int> > >::~SparseSolverBase()
  0.00     13.49     0.00        1     0.00     0.00  Eigen::MappedSparseMatrix<double, 0, int>::MappedSparseMatrix(long, long, long, int*, int*, double*, int*)
  0.00     13.49     0.00        1     0.00     0.00  Eigen::MappedSparseMatrix<double, 0, int>::~MappedSparseMatrix()
  0.00     13.49     0.00        1     0.00     0.00  Eigen::SparseCompressedBase<Eigen::Map<Eigen::SparseMatrix<double, 0, int>, 0, Eigen::Stride<0, 0> > >::SparseCompressedBase()
  0.00     13.49     0.00        1     0.00     0.00  Eigen::Map<Eigen::SparseMatrix<double, 0, int>, 0, Eigen::Stride<0, 0> >::Map(long, long, long, int*, int*, double*, int*)
  0.00     13.49     0.00        1     0.00     0.00  Eigen::Map<Eigen::SparseMatrix<double, 0, int>, 0, Eigen::Stride<0, 0> >::~Map()
  0.00     13.49     0.00        1     0.00     0.00  Eigen::Array<int, 2, 1, 0, 2, 1>::Array<int, int>(int const&, int const&)
  0.00     13.49     0.00        1     0.00     0.00  Eigen::SparseLU<Eigen::SparseMatrix<double, 0, int>, Eigen::COLAMDOrdering<int> >::initperfvalues()
  0.00     13.49     0.00        1     0.00     0.00  Eigen::SparseLU<Eigen::SparseMatrix<double, 0, int>, Eigen::COLAMDOrdering<int> >::SparseLU()
  0.00     13.49     0.00        1     0.00     0.00  Eigen::SparseLU<Eigen::SparseMatrix<double, 0, int>, Eigen::COLAMDOrdering<int> >::~SparseLU()
  0.00     13.49     0.00        1     0.00     0.00  Eigen::internal::plain_array<int, 2, 0, 0>::plain_array()
  0.00     13.49     0.00        1     0.00     0.00  Eigen::internal::LU_GlobalLU_t<Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<double, -1, 1, 0, -1, 1> >::LU_GlobalLU_t()
  0.00     13.49     0.00        1     0.00     0.00  Eigen::internal::LU_GlobalLU_t<Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<double, -1, 1, 0, -1, 1> >::~LU_GlobalLU_t()
  0.00     13.49     0.00        1     0.00     0.00  Eigen::internal::MappedSuperNodalMatrix<double, int>::MappedSuperNodalMatrix()
  0.00     13.49     0.00        1     0.00     0.00  Eigen::internal::MappedSuperNodalMatrix<double, int>::~MappedSuperNodalMatrix()
  0.00     13.49     0.00        1     0.00     0.00  void Eigen::internal::check_static_allocation_size<int, 2>()
  0.00     13.49     0.00        1     0.00     0.00  Eigen::ArrayBase<Eigen::Array<int, 2, 1, 0, 2, 1> >::ArrayBase()
  0.00     13.49     0.00        1     0.00     0.00  Eigen::DenseBase<Eigen::Array<int, 2, 1, 0, 2, 1> >::DenseBase()
  0.00     13.49     0.00        1     0.00     0.00  Cantera::ReactorBase::setInitialVolume(double)
  0.00     13.49     0.00        1     0.00     0.00  Cantera::MassFlowController::~MassFlowController()
  0.00     13.49     0.00        1     0.00     0.00  Cantera::PreconditionerBase::PreconditionerBase()
  0.00     13.49     0.00        1     0.00     0.00  Cantera::PreconditionerBase::~PreconditionerBase()
  0.00     13.49     0.00        1     0.00     0.00  Cantera::PressureController::setPressureCoeff(double)
  0.00     13.49     0.00        1     0.00     0.00  Cantera::PressureController::setMaster(Cantera::FlowDevice*)
  0.00     13.49     0.00        1     0.00     0.00  Cantera::PressureController::~PressureController()
  0.00     13.49     0.00        1     0.00     0.00  Cantera::ConstPressureReactor::ConstPressureReactor()
  0.00     13.49     0.00        1     0.00     0.00  Cantera::ConstPressureReactor::~ConstPressureReactor()
  0.00     13.49     0.00        1     0.00     0.00  Cantera::AdaptivePreconditioner::AdaptivePreconditioner()
  0.00     13.49     0.00        1     0.00     0.00  Cantera::AdaptivePreconditioner::~AdaptivePreconditioner()
  0.00     13.49     0.00        1     0.00     0.00  Cantera::IdealGasConstPressureReactor::IdealGasConstPressureReactor()
  0.00     13.49     0.00        1     0.00     0.00  Cantera::IdealGasConstPressureReactor::~IdealGasConstPressureReactor()
  0.00     13.49     0.00        1     0.00     0.00  Cantera::Reactor::~Reactor()
  0.00     13.49     0.00        1     0.00     2.50  Cantera::Solution::thermo()
  0.00     13.49     0.00        1     0.00     2.50  Cantera::Solution::kinetics()
  0.00     13.49     0.00        1     0.00     0.00  __gnu_cxx::new_allocator<Cantera::SensitivityParameter>::~new_allocator()
  0.00     13.49     0.00        1     0.00     0.00  __gnu_cxx::new_allocator<std::shared_ptr<Cantera::Solution> >::new_allocator()
  0.00     13.49     0.00        1     0.00     0.00  __gnu_cxx::new_allocator<std::shared_ptr<Cantera::Solution> >::~new_allocator()
  0.00     13.49     0.00        1     0.00     0.00  __gnu_cxx::new_allocator<double>::new_allocator()
  0.00     13.49     0.00        1     0.00     0.00  __gnu_cxx::new_allocator<unsigned long>::new_allocator()
  0.00     13.49     0.00        1     0.00     0.00  std::__shared_ptr<Cantera::Kinetics, (__gnu_cxx::_Lock_policy)2>::get() const
  0.00     13.49     0.00        1     0.00     0.00  std::__shared_ptr_access<Cantera::ThermoPhase, (__gnu_cxx::_Lock_policy)2, false, false>::operator->() const
  0.00     13.49     0.00        1     0.00     0.00  std::__shared_ptr_access<Cantera::Kinetics, (__gnu_cxx::_Lock_policy)2, false, false>::_M_get() const
  0.00     13.49     0.00        1     0.00     0.00  std::__shared_ptr_access<Cantera::Kinetics, (__gnu_cxx::_Lock_policy)2, false, false>::operator*() const
  0.00     13.49     0.00        1     0.00     0.00  std::allocator<Cantera::SensitivityParameter>::~allocator()
  0.00     13.49     0.00        1     0.00     0.00  std::allocator<std::shared_ptr<Cantera::Solution> >::allocator()
  0.00     13.49     0.00        1     0.00     0.00  std::allocator<std::shared_ptr<Cantera::Solution> >::~allocator()
  0.00     13.49     0.00        1     0.00     0.00  std::allocator<double>::allocator()
  0.00     13.49     0.00        1     0.00     0.00  std::allocator<unsigned long>::allocator()
  0.00     13.49     0.00        1     0.00     2.50  std::shared_ptr<Cantera::ThermoPhase>::shared_ptr(std::shared_ptr<Cantera::ThermoPhase> const&)
  0.00     13.49     0.00        1     0.00     0.00  std::shared_ptr<Cantera::ThermoPhase>::~shared_ptr()
  0.00     13.49     0.00        1     0.00     2.50  std::shared_ptr<Cantera::Kinetics>::shared_ptr(std::shared_ptr<Cantera::Kinetics> const&)
  0.00     13.49     0.00        1     0.00     0.00  std::shared_ptr<Cantera::Kinetics>::~shared_ptr()
  0.00     13.49     0.00        1     0.00     0.00  std::shared_ptr<Cantera::Solution>::~shared_ptr()
  0.00     13.49     0.00        1     0.00     0.00  void std::_Destroy_aux<false>::__destroy<std::shared_ptr<Cantera::Solution>*>(std::shared_ptr<Cantera::Solution>*, std::shared_ptr<Cantera::Solution>*)
  0.00     13.49     0.00        1     0.00     0.00  void std::_Destroy_aux<true>::__destroy<Cantera::SensitivityParameter*>(Cantera::SensitivityParameter*, Cantera::SensitivityParameter*)
  0.00     13.49     0.00        1     0.00     0.00  std::_Vector_base<Cantera::SensitivityParameter, std::allocator<Cantera::SensitivityParameter> >::_Vector_impl::~_Vector_impl()
  0.00     13.49     0.00        1     0.00     0.00  std::_Vector_base<Cantera::SensitivityParameter, std::allocator<Cantera::SensitivityParameter> >::_M_deallocate(Cantera::SensitivityParameter*, unsigned long)
  0.00     13.49     0.00        1     0.00     0.00  std::_Vector_base<Cantera::SensitivityParameter, std::allocator<Cantera::SensitivityParameter> >::_M_get_Tp_allocator()
  0.00     13.49     0.00        1     0.00     0.00  std::_Vector_base<Cantera::SensitivityParameter, std::allocator<Cantera::SensitivityParameter> >::~_Vector_base()
  0.00     13.49     0.00        1     0.00     0.00  std::_Vector_base<std::shared_ptr<Cantera::Solution>, std::allocator<std::shared_ptr<Cantera::Solution> > >::_Vector_impl::_Vector_impl()
  0.00     13.49     0.00        1     0.00     0.00  std::_Vector_base<std::shared_ptr<Cantera::Solution>, std::allocator<std::shared_ptr<Cantera::Solution> > >::_Vector_impl::~_Vector_impl()
  0.00     13.49     0.00        1     0.00     0.00  std::_Vector_base<std::shared_ptr<Cantera::Solution>, std::allocator<std::shared_ptr<Cantera::Solution> > >::_M_deallocate(std::shared_ptr<Cantera::Solution>*, unsigned long)
  0.00     13.49     0.00        1     0.00     0.00  std::_Vector_base<std::shared_ptr<Cantera::Solution>, std::allocator<std::shared_ptr<Cantera::Solution> > >::_Vector_impl_data::_Vector_impl_data()
  0.00     13.49     0.00        1     0.00     0.00  std::_Vector_base<std::shared_ptr<Cantera::Solution>, std::allocator<std::shared_ptr<Cantera::Solution> > >::_M_get_Tp_allocator()
  0.00     13.49     0.00        1     0.00     0.00  std::_Vector_base<std::shared_ptr<Cantera::Solution>, std::allocator<std::shared_ptr<Cantera::Solution> > >::_Vector_base()
  0.00     13.49     0.00        1     0.00     0.00  std::_Vector_base<std::shared_ptr<Cantera::Solution>, std::allocator<std::shared_ptr<Cantera::Solution> > >::~_Vector_base()
  0.00     13.49     0.00        1     0.00     0.00  std::_Vector_base<double, std::allocator<double> >::_Vector_impl::_Vector_impl()
  0.00     13.49     0.00        1     0.00     0.00  std::_Vector_base<double, std::allocator<double> >::_Vector_impl_data::_Vector_impl_data()
  0.00     13.49     0.00        1     0.00     0.00  std::_Vector_base<double, std::allocator<double> >::_Vector_base()
  0.00     13.49     0.00        1     0.00     0.00  std::_Vector_base<unsigned long, std::allocator<unsigned long> >::_Vector_impl::_Vector_impl()
  0.00     13.49     0.00        1     0.00     0.00  std::_Vector_base<unsigned long, std::allocator<unsigned long> >::_Vector_impl_data::_Vector_impl_data()
  0.00     13.49     0.00        1     0.00     0.00  std::_Vector_base<unsigned long, std::allocator<unsigned long> >::_Vector_base()
  0.00     13.49     0.00        1     0.00     2.50  std::__shared_ptr<Cantera::ThermoPhase, (__gnu_cxx::_Lock_policy)2>::__shared_ptr(std::__shared_ptr<Cantera::ThermoPhase, (__gnu_cxx::_Lock_policy)2> const&)
  0.00     13.49     0.00        1     0.00     0.00  std::__shared_ptr<Cantera::ThermoPhase, (__gnu_cxx::_Lock_policy)2>::~__shared_ptr()
  0.00     13.49     0.00        1     0.00     2.50  std::__shared_ptr<Cantera::Kinetics, (__gnu_cxx::_Lock_policy)2>::__shared_ptr(std::__shared_ptr<Cantera::Kinetics, (__gnu_cxx::_Lock_policy)2> const&)
  0.00     13.49     0.00        1     0.00     0.00  std::__shared_ptr<Cantera::Kinetics, (__gnu_cxx::_Lock_policy)2>::~__shared_ptr()
  0.00     13.49     0.00        1     0.00     0.00  std::__shared_ptr<Cantera::Solution, (__gnu_cxx::_Lock_policy)2>::~__shared_ptr()
  0.00     13.49     0.00        1     0.00     0.00  std::vector<Cantera::SensitivityParameter, std::allocator<Cantera::SensitivityParameter> >::~vector()
  0.00     13.49     0.00        1     0.00     0.00  std::vector<std::shared_ptr<Cantera::Solution>, std::allocator<std::shared_ptr<Cantera::Solution> > >::vector()
  0.00     13.49     0.00        1     0.00     0.00  std::vector<std::shared_ptr<Cantera::Solution>, std::allocator<std::shared_ptr<Cantera::Solution> > >::~vector()
  0.00     13.49     0.00        1     0.00     0.00  std::vector<double, std::allocator<double> >::vector()
  0.00     13.49     0.00        1     0.00     0.00  std::vector<unsigned long, std::allocator<unsigned long> >::vector()
  0.00     13.49     0.00        1     0.00     0.00  void std::_Destroy<Cantera::SensitivityParameter*>(Cantera::SensitivityParameter*, Cantera::SensitivityParameter*)
  0.00     13.49     0.00        1     0.00     0.00  void std::_Destroy<Cantera::SensitivityParameter*, Cantera::SensitivityParameter>(Cantera::SensitivityParameter*, Cantera::SensitivityParameter*, std::allocator<Cantera::SensitivityParameter>&)
  0.00     13.49     0.00        1     0.00     0.00  void std::_Destroy<std::shared_ptr<Cantera::Solution>*>(std::shared_ptr<Cantera::Solution>*, std::shared_ptr<Cantera::Solution>*)
  0.00     13.49     0.00        1     0.00     0.00  void std::_Destroy<std::shared_ptr<Cantera::Solution>*, std::shared_ptr<Cantera::Solution> >(std::shared_ptr<Cantera::Solution>*, std::shared_ptr<Cantera::Solution>*, std::allocator<std::shared_ptr<Cantera::Solution> >&)
