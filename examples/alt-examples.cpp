class ReactType;
class IdealReactType;
class ConditionerBase
  {
  public:
    ConditionerBase(/* args */){};
    ~ConditionerBase(){};
    virtual void visit(ReactType* reactor);
    virtual void visit(IdealReactType* reactor);
  };

  class ReactType
  {
  public:
    ReactType(/* args */){};
    ~ReactType(){};
    virtual void accept(ConditionerBase *conditioner);
  };
  
  class IdealReactType: public ReactType
  {
  public:
    IdealReactType(/* args */){};
    ~IdealReactType(){};
    virtual void accept(ConditionerBase *conditioner);
  };

  class Conditioner : public ConditionerBase
  {
  public:
    Conditioner(/* args */){};
    ~Conditioner(){};
    virtual void visit(ReactType* reactor){std::cout<<"Visited unknown reactor"<<std::endl;};
    virtual void visit(IdealReactType* ideal_reactor){std::cout<<"Visited ideal_reactor"<<std::endl;};
  };

  void ReactType::accept(ConditionerBase *conditioner){conditioner->visit(this);};
  void IdealReactType::accept(ConditionerBase *conditioner){conditioner->visit(this);};
  void ConditionerBase::visit(ReactType* reactor){std::cout<<"Base Reactor"<<std::endl;};
  void ConditionerBase::visit(IdealReactType* reactor){std::cout<<"Base Ideal"<<std::endl;};

void doubleDispatchExample()
{

  std::vector<ReactType*> reactors;
  ReactType* unknown = new ReactType();
  IdealReactType* ideal_reactor = new IdealReactType();
  reactors.push_back(unknown);
  reactors.push_back(ideal_reactor);
  Conditioner conditioner;
  for (size_t i = 0; i < 2; i++)
  {
    reactors.at(i)->accept(&conditioner);
  }
  delete unknown;
  delete ideal_reactor;
}

void memcpyExample()
{
  int size = 4;
  double *b  = new double[size];
  double *x  = new double[size];

  for (size_t i = 0; i < size; i++)
  {
    b[i] = i;
  }
  

  std::memcpy(x,b,size*sizeof(b));

  for (size_t i = 0; i < size; i++)
  {
    std::cout<<x[i]<<std::endl;
  }
}