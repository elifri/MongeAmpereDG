/// taken from http://www.breaktrycatch.com/deeper-into-the-rabbit-hole-c-function-pointers-and-callbacks/

/*********************************
*Class: Function
*Description: Based off of Elbert Mai's (http://www.codeproject.com/script/Membership/View.aspx?mid=2301380) implementation
*of Callbacks (http://www.codeproject.com/KB/cpp/CPPCallback.aspx).
*Modified to support equality and heavily commented to provide easier understanding of what's happening under the hood.
*Author: jkeon
**********************************/
 
#ifndef _FUNCTION_H_
#define _FUNCTION_H_
 
//////////////////////////////////////////////////////////////////////
// MACROS ////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
 
//Creates and Returns a Free Function object. Static functions etc.
#define FREE_FUNCTION(functionPointer) util::CreateFunctionBuilder(functionPointer).Wrap<functionPointer>()
 
//Creates and Returns a Member Function object. Functions inside a class.
#define MEMBER_FUNCTION(functionPointer, instancePointer) util::CreateFunctionBuilder(functionPointer).Wrap<functionPointer>(instancePointer)
 
//////////////////////////////////////////////////////////////////////
// NAMESPACE /////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
 
namespace util {
 
    //////////////////////////////////////////////////////////////////////
    // CLASS DECLARATION /////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////
 
    //Declaring a Class called Function which can be of any type. Naturally we want it to be the Function Signature.
    //We're declaring this here so that the Type itself exists but does nothing and we can extend this with Partial Template Specialization later.
    template <typename FunctionSignature>
    class Function;
 
    //////////////////////////////////////////////////////////////////////
    // GLOBAL EQUALITY CHECKS ////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////
 
    template <typename EqualFuncReturnType>
    inline bool operator== (const Function<EqualFuncReturnType> &lhs, const Function<EqualFuncReturnType> &rhs) {
        return ((lhs.functionPointer == rhs.functionPointer) && (lhs.instancePointer == rhs.instancePointer));
    }
 
    template <typename NotEqualFuncReturnType>
    inline bool operator!= (const Function<NotEqualFuncReturnType> &lhs, const Function<NotEqualFuncReturnType> &rhs) {
        return ((lhs.functionPointer != rhs.functionPointer) || (lhs.instancePointer != rhs.instancePointer));
    }
 
    //*********************************************************************************************************************************************************************
    //*********************************************************************************************************************************************************************
 
    //////////////////////////////////////////////////////////////////////
    // FUNCTION: 0 PARAMETER VERSION /////////////////////////////////////
    //////////////////////////////////////////////////////////////////////
 
    //Using Template Partial Specialization we will only allow implementations of a Function class which are typed to be a Function that could have any return type.
    template <typename ReturnType>
    class Function<ReturnType ()> {
 
    //PUBLIC FUNCTIONS
    public:
        //Default Constructor - To allow for creation with assigning
        inline Function() : functionPointer(0), instancePointer(0) {}
 
        //Constructor - Set function pointer via initalizer list
        inline Function(ReturnType (*fp)(const void*), const void *ip) : functionPointer(fp), instancePointer(ip) {}
 
        //Copy Constructor - Set function pointer via initializer list
        inline Function(const Function& other) : functionPointer(other.functionPointer), instancePointer(other.instancePointer) {}
 
        //Destructor - Simply set function and instance pointer to 0
        inline ~Function() {
            functionPointer = 0;
            instancePointer = 0;
        }
 
        //Overloading = operator to allow for assignment
        inline Function<ReturnType>& operator= (const Function<ReturnType> &other) {
            functionPointer = other.functionPointer;
            instancePointer = other.instancePointer;
            return *this;
        }
 
        //Overloading () operator so we can use this like a Function
        inline ReturnType operator() () const {
            return (*functionPointer)(instancePointer);
        }
 
        //Safe Bool Idiom - Allows for checking if a Function has a valid functionPointer or not.
        typedef const void* Function::*bool_type;
        inline operator bool_type() const {
            return (functionPointer != 0) ? &Function::instancePointer : false;
        }
 
        //Allowing for Equality/InEquality Checks by granting the Global Equality/InEquality Check Functions access to this classes internals via the friend keyword.
        template <typename EqualFuncReturnType>
        friend inline bool operator ==(const Function<EqualFuncReturnType> & lhs, const Function<EqualFuncReturnType> & rhs);
        template<typename NotEqualFuncReturnType>
        friend inline bool operator !=(const Function<NotEqualFuncReturnType> & lhs, const Function<NotEqualFuncReturnType> & rhs);
 
    //PRIVATE VARIABLES
    private:
        //Stores a Free Function Pointer to the static wrapper function
        ReturnType (*functionPointer)(const void*);
        //Stores an Instance Pointer
        const void *instancePointer;
 
    };
 
    //////////////////////////////////////////////////////////////////////
    // FUNCTION BUILDER: FREE FUNCTION 0 PARAMETER VERSION ///////////////
    //////////////////////////////////////////////////////////////////////
 
    //Builds a Free Function
    template <typename ReturnType>
    class FreeFunctionBuilder0 {
    public:
        //Empty Constructor/Destructor
        inline FreeFunctionBuilder0() {}
        inline ~FreeFunctionBuilder0() {}
 
        //Performs the wrapping from the actual Free Function to the static Free Function Wrapper in this class.
        template<ReturnType (*functionPointer)()>
        inline static Function<ReturnType ()> Wrap() {
            return Function<ReturnType ()>(&FreeFunctionBuilder0::template Wrapper<functionPointer>, 0);
        }
 
    private:
        //Redirects to the functionPointer passed in at compile time via template params
        template<ReturnType (*functionPointer)()>
        inline static ReturnType Wrapper(const void*) {
            return (*functionPointer)();
        }
    };
 
    //////////////////////////////////////////////////////////////////////
    // FUNCTION BUILDER: MEMBER FUNCTION 0 PARAMETER VERSION /////////////
    //////////////////////////////////////////////////////////////////////
 
    //Builds a Member Function
    template <typename ReturnType, class ClassType>
    class MemberFunctionBuilder0 {
    public:
        //Empty Constructor/Destructor
        inline MemberFunctionBuilder0() {}
        inline ~MemberFunctionBuilder0() {}
 
        //Performs the wrapping from the actual Member Function to the static Free Function Wrapper in this class. Casts the instance pointer to a const void*
        template<ReturnType (ClassType::*functionPointer)()>
        inline static Function<ReturnType ()> Wrap(ClassType* ip) {
            return Function<ReturnType ()>(&MemberFunctionBuilder0::template Wrapper<functionPointer>, static_cast<const void*>(ip));
        }
 
    private:
        //Redirects to the functionPointer passed in at compile time via template params. Also casts the const void* back to the correct class instance.
        template<ReturnType (ClassType::*functionPointer)()>
        inline static ReturnType Wrapper(const void* ip) {
            ClassType* instancePointer = const_cast<ClassType*>(static_cast<const ClassType*>(ip));
            return (instancePointer->*functionPointer)();
        }
 
    };
 
    //////////////////////////////////////////////////////////////////////
    // FUNCTION BUILDER: CONST MEMBER FUNCTION 0 PARAMETER VERSION ///////
    //////////////////////////////////////////////////////////////////////
 
    //Builds a Const Member Function
    template <typename ReturnType, class ClassType>
    class ConstMemberFunctionBuilder0 {
    public:
        //Empty Constructor/Destructor
        inline ConstMemberFunctionBuilder0() {}
        inline ~ConstMemberFunctionBuilder0() {}
 
        //Performs the wrapping from the actual Const Member Function to the static Free Function Wrapper in this class. Casts the instance pointer to a const void*
        template<ReturnType (ClassType::*functionPointer)() const>
        inline static Function<ReturnType ()> Wrap(ClassType* ip) {
            return Function<ReturnType ()>(&ConstMemberFunctionBuilder0::template Wrapper<functionPointer>, static_cast<const void*>(ip));
        }
 
    private:
        //Redirects to the functionPointer passed in at compile time via template params. Also casts the const void* back to the correct class instance.
        template<ReturnType (ClassType::*functionPointer)() const>
        inline static ReturnType Wrapper(const void* ip) {
            ClassType* instancePointer = const_cast<ClassType*>(static_cast<const ClassType*>(ip));
            return (instancePointer->*functionPointer)();
        }
 
    };
 
    //////////////////////////////////////////////////////////////////////
    // GLOBAL FUNCTION CREATOR: FREE FUNCTION 0 PARAMETER VERSION ////////
    //////////////////////////////////////////////////////////////////////
 
    template <typename ReturnType>
    inline FreeFunctionBuilder0<ReturnType> CreateFunctionBuilder(ReturnType (*fp)()) {
        return FreeFunctionBuilder0<ReturnType>();
    }
 
    //////////////////////////////////////////////////////////////////////
    // GLOBAL FUNCTION CREATOR: MEMBER FUNCTION 0 PARAMETER VERSION //////
    //////////////////////////////////////////////////////////////////////
 
    template <typename ReturnType, class ClassType>
    inline MemberFunctionBuilder0<ReturnType, ClassType> CreateFunctionBuilder(ReturnType (ClassType::*fp)()) {
        return MemberFunctionBuilder0<ReturnType, ClassType>();
    }
 
    //////////////////////////////////////////////////////////////////////
    // GLOBAL FUNCTION CREATOR: C MEMBER FUNCTION 0 PARAMETER VERSION ////
    //////////////////////////////////////////////////////////////////////
 
    template <typename ReturnType, class ClassType>
    inline ConstMemberFunctionBuilder0<ReturnType, ClassType> CreateFunctionBuilder(ReturnType (ClassType::*fp)() const) {
        return ConstMemberFunctionBuilder0<ReturnType, ClassType>();
    }
 
    //*********************************************************************************************************************************************************************
    //*********************************************************************************************************************************************************************
 
    //////////////////////////////////////////////////////////////////////
    // FUNCTION: 1 PARAMETER VERSION /////////////////////////////////////
    //////////////////////////////////////////////////////////////////////
 
    //Using Template Partial Specialization we will only allow implementations of a Function class which are typed to be a Function that could have any return type.
    template <typename ReturnType, typename Param0>
    class Function<ReturnType (Param0)> {
 
    //PUBLIC FUNCTIONS
    public:
        //Default Constructor - To allow for creation with assigning
        inline Function() : functionPointer(0), instancePointer(0) {}
 
        //Constructor - Set function pointer via initalizer list
        inline Function(ReturnType (*fp)(const void*, Param0 p0), const void *ip) : functionPointer(fp), instancePointer(ip) {}
 
        //Copy Constructor - Set function pointer via initializer list
        inline Function(const Function& other) : functionPointer(other.functionPointer), instancePointer(other.instancePointer) {}
 
        //Destructor - Simply set function and instance pointer to 0
        inline ~Function() {
            functionPointer = 0;
            instancePointer = 0;
        }
 
        //Overloading = operator to allow for assignment
        inline Function<ReturnType>& operator= (const Function<ReturnType> &other) {
            functionPointer = other.functionPointer;
            instancePointer = other.instancePointer;
            return *this;
        }
 
        //Overloading () operator so we can use this like a Function
        inline ReturnType operator() (Param0 p0) const {
            return (*functionPointer)(instancePointer, p0);
        }
 
        //Safe Bool Idiom - Allows for checking if a Function has a valid functionPointer or not.
        typedef const void* Function::*bool_type;
        inline operator bool_type() const {
            return (functionPointer != 0) ? &Function::instancePointer : false;
        }
 
        //Allowing for Equality/InEquality Checks by granting the Global Equality/InEquality Check Functions access to this classes internals via the friend keyword.
        template <typename EqualFuncReturnType>
        friend inline bool operator ==(const Function<EqualFuncReturnType> & lhs, const Function<EqualFuncReturnType> & rhs);
        template<typename NotEqualFuncReturnType>
        friend inline bool operator !=(const Function<NotEqualFuncReturnType> & lhs, const Function<NotEqualFuncReturnType> & rhs);
 
    //PRIVATE VARIABLES
    private:
        //Stores a Free Function Pointer to the static wrapper function
        ReturnType (*functionPointer)(const void*, Param0 p0);
        //Stores an Instance Pointer
        const void *instancePointer;
 
    };
 
    //////////////////////////////////////////////////////////////////////
    // FUNCTION BUILDER: FREE FUNCTION 1 PARAMETER VERSION ///////////////
    //////////////////////////////////////////////////////////////////////
 
    //Builds a Free Function
    template <typename ReturnType, typename Param0>
    class FreeFunctionBuilder1 {
    public:
        //Empty Constructor/Destructor
        inline FreeFunctionBuilder1() {}
        inline ~FreeFunctionBuilder1() {}
 
        //Performs the wrapping from the actual Free Function to the static Free Function Wrapper in this class.
        template<ReturnType (*functionPointer)(Param0 p0)>
        inline static Function<ReturnType (Param0)> Wrap() {
            return Function<ReturnType (Param0)>(&FreeFunctionBuilder1::template Wrapper<functionPointer>, 0);
        }
 
    private:
        //Redirects to the functionPointer passed in at compile time via template params
        template<ReturnType (*functionPointer)(Param0 p0)>
        inline static ReturnType Wrapper(const void*, Param0 p0) {
            return (*functionPointer)(p0);
        }
    };
 
    //////////////////////////////////////////////////////////////////////
    // FUNCTION BUILDER: MEMBER FUNCTION 1 PARAMETER VERSION /////////////
    //////////////////////////////////////////////////////////////////////
 
    //Builds a Member Function
    template <typename ReturnType, class ClassType, typename Param0>
    class MemberFunctionBuilder1 {
    public:
        //Empty Constructor/Destructor
        inline MemberFunctionBuilder1() {}
        inline ~MemberFunctionBuilder1() {}
 
        //Performs the wrapping from the actual Member Function to the static Free Function Wrapper in this class. Casts the instance pointer to a const void*
        template<ReturnType (ClassType::*functionPointer)(Param0 p0)>
        inline static Function<ReturnType (Param0)> Wrap(ClassType* ip) {
            return Function<ReturnType (Param0)>(&MemberFunctionBuilder1::template Wrapper<functionPointer>, static_cast<const void*>(ip));
        }
 
    private:
        //Redirects to the functionPointer passed in at compile time via template params. Also casts the const void* back to the correct class instance.
        template<ReturnType (ClassType::*functionPointer)(Param0 p0)>
        inline static ReturnType Wrapper(const void* ip, Param0 p0) {
            ClassType* instancePointer = const_cast<ClassType*>(static_cast<const ClassType*>(ip));
            return (instancePointer->*functionPointer)(p0);
        }
 
    };
 
    //////////////////////////////////////////////////////////////////////
    // FUNCTION BUILDER: CONST MEMBER FUNCTION 1 PARAMETER VERSION ///////
    //////////////////////////////////////////////////////////////////////
 
    //Builds a Const Member Function
    template <typename ReturnType, class ClassType, typename Param0>
    class ConstMemberFunctionBuilder1 {
    public:
        //Empty Constructor/Destructor
        inline ConstMemberFunctionBuilder1() {}
        inline ~ConstMemberFunctionBuilder1() {}
 
        //Performs the wrapping from the actual Const Member Function to the static Free Function Wrapper in this class. Casts the instance pointer to a const void*
        template<ReturnType (ClassType::*functionPointer)(Param0 p0) const>
        inline static Function<ReturnType (Param0 p0)> Wrap(ClassType* ip) {
            return Function<ReturnType (Param0 p0)>(&ConstMemberFunctionBuilder1::template Wrapper<functionPointer>, static_cast<const void*>(ip));
        }
 
    private:
        //Redirects to the functionPointer passed in at compile time via template params. Also casts the const void* back to the correct class instance.
        template<ReturnType (ClassType::*functionPointer)(Param0 p0) const>
        inline static ReturnType Wrapper(const void* ip, Param0 p0) {
            ClassType* instancePointer = const_cast<ClassType*>(static_cast<const ClassType*>(ip));
            return (instancePointer->*functionPointer)(p0);
        }
 
    };
 
    //////////////////////////////////////////////////////////////////////
    // GLOBAL FUNCTION CREATOR: FREE FUNCTION 1 PARAMETER VERSION ////////
    //////////////////////////////////////////////////////////////////////
 
    template <typename ReturnType, typename Param0>
    inline FreeFunctionBuilder1<ReturnType, Param0> CreateFunctionBuilder(ReturnType (*fp)(Param0 p0)) {
        return FreeFunctionBuilder1<ReturnType, Param0>();
    }
 
    //////////////////////////////////////////////////////////////////////
    // GLOBAL FUNCTION CREATOR: MEMBER FUNCTION 1 PARAMETER VERSION //////
    //////////////////////////////////////////////////////////////////////
 
    template <typename ReturnType, class ClassType, typename Param0>
    inline MemberFunctionBuilder1<ReturnType, ClassType, Param0> CreateFunctionBuilder(ReturnType (ClassType::*fp)(Param0 p0)) {
        return MemberFunctionBuilder1<ReturnType, ClassType, Param0>();
    }
 
    //////////////////////////////////////////////////////////////////////
    // GLOBAL FUNCTION CREATOR: C MEMBER FUNCTION 1 PARAMETER VERSION ////
    //////////////////////////////////////////////////////////////////////
 
    template <typename ReturnType, class ClassType, typename Param0>
    inline ConstMemberFunctionBuilder1<ReturnType, ClassType, Param0> CreateFunctionBuilder(ReturnType (ClassType::*fp)(Param0 p0) const) {
        return ConstMemberFunctionBuilder1<ReturnType, ClassType, Param0>();
    }
 
    //*********************************************************************************************************************************************************************
    //*********************************************************************************************************************************************************************
 
    //////////////////////////////////////////////////////////////////////
    // FUNCTION: 2 PARAMETER VERSION /////////////////////////////////////
    //////////////////////////////////////////////////////////////////////
 
    //Using Template Partial Specialization we will only allow implementations of a Function class which are typed to be a Function that could have any return type.
    template <typename ReturnType, typename Param0, typename Param1>
    class Function<ReturnType (Param0, Param1)> {
 
    //PUBLIC FUNCTIONS
    public:
        //Default Constructor - To allow for creation with assigning
        inline Function() : functionPointer(0), instancePointer(0) {}
 
        //Constructor - Set function pointer via initalizer list
        inline Function(ReturnType (*fp)(const void*, Param0 p0, Param1 p1), const void *ip) : functionPointer(fp), instancePointer(ip) {}
 
        //Copy Constructor - Set function pointer via initializer list
        inline Function(const Function& other) : functionPointer(other.functionPointer), instancePointer(other.instancePointer) {}
 
        //Destructor - Simply set function and instance pointer to 0
        inline ~Function() {
            functionPointer = 0;
            instancePointer = 0;
        }
 
        //Overloading = operator to allow for assignment
        inline Function<ReturnType>& operator= (const Function<ReturnType> &other) {
            functionPointer = other.functionPointer;
            instancePointer = other.instancePointer;
            return *this;
        }
 
        //Overloading () operator so we can use this like a Function
        inline ReturnType operator() (Param0 p0, Param1 p1) const {
            return (*functionPointer)(instancePointer, p0, p1);
        }
 
        //Safe Bool Idiom - Allows for checking if a Function has a valid functionPointer or not.
        typedef const void* Function::*bool_type;
        inline operator bool_type() const {
            return (functionPointer != 0) ? &Function::instancePointer : false;
        }
 
        template <typename EqualFuncReturnType>
        friend inline bool operator ==(const Function<EqualFuncReturnType> & lhs, const Function<EqualFuncReturnType> & rhs);
        template<typename NotEqualFuncReturnType>
        friend inline bool operator !=(const Function<NotEqualFuncReturnType> & lhs, const Function<NotEqualFuncReturnType> & rhs);
 
    //PRIVATE VARIABLES
    private:
        //Stores a Free Function Pointer to the static wrapper function
        ReturnType (*functionPointer)(const void*, Param0 p0, Param1 p1);
        //Stores an Instance Pointer
        const void *instancePointer;
 
    };
 
    //////////////////////////////////////////////////////////////////////
    // FUNCTION BUILDER: FREE FUNCTION 2 PARAMETER VERSION ///////////////
    //////////////////////////////////////////////////////////////////////
 
    //Builds a Free Function
    template <typename ReturnType, typename Param0, typename Param1>
    class FreeFunctionBuilder2 {
    public:
        //Empty Constructor/Destructor
        inline FreeFunctionBuilder2() {}
        inline ~FreeFunctionBuilder2() {}
 
        //Performs the wrapping from the actual Free Function to the static Free Function Wrapper in this class.
        template<ReturnType (*functionPointer)(Param0 p0, Param1 p1)>
        inline static Function<ReturnType (Param0, Param1)> Wrap() {
            return Function<ReturnType (Param0, Param1)>(&FreeFunctionBuilder2::template Wrapper<functionPointer>, 0);
        }
 
    private:
        //Redirects to the functionPointer passed in at compile time via template params
        template<ReturnType (*functionPointer)(Param0 p0, Param1 p1)>
        inline static ReturnType Wrapper(const void*, Param0 p0, Param1 p1) {
            return (*functionPointer)(p0, p1);
        }
    };
 
    //////////////////////////////////////////////////////////////////////
    // FUNCTION BUILDER: MEMBER FUNCTION 2 PARAMETER VERSION /////////////
    //////////////////////////////////////////////////////////////////////
 
    //Builds a Member Function
    template <typename ReturnType, class ClassType, typename Param0, typename Param1>
    class MemberFunctionBuilder2 {
    public:
        //Empty Constructor/Destructor
        inline MemberFunctionBuilder2() {}
        inline ~MemberFunctionBuilder2() {}
 
        //Performs the wrapping from the actual Member Function to the static Free Function Wrapper in this class. Casts the instance pointer to a const void*
        template<ReturnType (ClassType::*functionPointer)(Param0 p0, Param1 p1)>
        inline static Function<ReturnType (Param0, Param1)> Wrap(ClassType* ip) {
            return Function<ReturnType (Param0, Param1)>(&MemberFunctionBuilder2::template Wrapper<functionPointer>, static_cast<const void*>(ip));
        }
 
    private:
        //Redirects to the functionPointer passed in at compile time via template params. Also casts the const void* back to the correct class instance.
        template<ReturnType (ClassType::*functionPointer)(Param0 p0, Param1 p1)>
        inline static ReturnType Wrapper(const void* ip, Param0 p0, Param1 p1) {
            ClassType* instancePointer = const_cast<ClassType*>(static_cast<const ClassType*>(ip));
            return (instancePointer->*functionPointer)(p0, p1);
        }
 
    };
 
    //////////////////////////////////////////////////////////////////////
    // FUNCTION BUILDER: CONST MEMBER FUNCTION 2 PARAMETER VERSION ///////
    //////////////////////////////////////////////////////////////////////
 
    //Builds a Const Member Function
    template <typename ReturnType, class ClassType, typename Param0, typename Param1>
    class ConstMemberFunctionBuilder2 {
    public:
        //Empty Constructor/Destructor
        inline ConstMemberFunctionBuilder2() {}
        inline ~ConstMemberFunctionBuilder2() {}
 
        //Performs the wrapping from the actual Const Member Function to the static Free Function Wrapper in this class. Casts the instance pointer to a const void*
        template<ReturnType (ClassType::*functionPointer)(Param0 p0, Param0 p1) const>
        inline static Function<ReturnType (Param0 p0, Param1 p1)> Wrap(ClassType* ip) {
            return Function<ReturnType (Param0 p0, Param1 p1)>(&ConstMemberFunctionBuilder2::template Wrapper<functionPointer>, static_cast<const void*>(ip));
        }
 
    private:
        //Redirects to the functionPointer passed in at compile time via template params. Also casts the const void* back to the correct class instance.
        template<ReturnType (ClassType::*functionPointer)(Param0 p0, Param1 p1) const>
        inline static ReturnType Wrapper(const void* ip, Param0 p0, Param1 p1) {
            ClassType* instancePointer = const_cast<ClassType*>(static_cast<const ClassType*>(ip));
            return (instancePointer->*functionPointer)(p0, p1);
        }
 
    };
 
    //////////////////////////////////////////////////////////////////////
    // GLOBAL FUNCTION CREATOR: FREE FUNCTION 2 PARAMETER VERSION ////////
    //////////////////////////////////////////////////////////////////////
 
    template <typename ReturnType, typename Param0, typename Param1>
    inline FreeFunctionBuilder2<ReturnType, Param0, Param1> CreateFunctionBuilder(ReturnType (*fp)(Param0 p0, Param1 p1)) {
        return FreeFunctionBuilder2<ReturnType, Param0, Param1>();
    }
 
    //////////////////////////////////////////////////////////////////////
    // GLOBAL FUNCTION CREATOR: MEMBER FUNCTION 2 PARAMETER VERSION //////
    //////////////////////////////////////////////////////////////////////
 
    template <typename ReturnType, class ClassType, typename Param0, typename Param1>
    inline MemberFunctionBuilder2<ReturnType, ClassType, Param0, Param1> CreateFunctionBuilder(ReturnType (ClassType::*fp)(Param0 p0, Param1 p1)) {
        return MemberFunctionBuilder2<ReturnType, ClassType, Param0, Param1>();
    }
 
    //////////////////////////////////////////////////////////////////////
    // GLOBAL FUNCTION CREATOR: C MEMBER FUNCTION 2 PARAMETER VERSION ////
    //////////////////////////////////////////////////////////////////////
 
    template <typename ReturnType, class ClassType, typename Param0, typename Param1>
    inline ConstMemberFunctionBuilder2<ReturnType, ClassType, Param0, Param1> CreateFunctionBuilder(ReturnType (ClassType::*fp)(Param0 p0, Param1 p1) const) {
        return ConstMemberFunctionBuilder2<ReturnType, ClassType, Param0, Param1>();
    }
}
#endif
