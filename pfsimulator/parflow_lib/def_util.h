#ifndef _DEF_UTIL_H
#define _DEF_UTIL_H
/*--------------------------------------------------------------------------------
  Header file containing all the macro utilities needed.
  Here be dragons.
  --------------------------------------------------------------------------------*/

/* Defer pragma insertions, use in nested macro chains */
#define PRAGMA(args) _Pragma(args)

/* Evaluator macros, use to evaluate macro chains */
#define EVAL(...) EVAL1024(__VA_ARGS__)
#define EVAL1024(...) EVAL512(EVAL512(__VA_ARGS__))
#define EVAL512(...) EVAL256(EVAL256(__VA_ARGS__))
#define EVAL256(...) EVAL128(EVAL128(__VA_ARGS__))
#define EVAL128(...) EVAL64(EVAL64(__VA_ARGS__))
#define EVAL64(...) EVAL32(EVAL32(__VA_ARGS__))
#define EVAL32(...) EVAL16(EVAL16(__VA_ARGS__))
#define EVAL16(...) EVAL8(EVAL8(__VA_ARGS__))
#define EVAL8(...) EVAL4(EVAL4(__VA_ARGS__))
#define EVAL4(...) EVAL2(EVAL2(__VA_ARGS__))
#define EVAL2(...) EVAL1(EVAL1(__VA_ARGS__))
#define EVAL1(...) __VA_ARGS__

/* Does what it says on the tin */
#define FIRST(a, ...) a
#define SECOND(a, b, ...) b
#define _SECOND(...) DEFER(SECOND)(__VA_ARGS__)
#define REST(a, b, ...) __VA_ARGS__
#define _REST(...) DEFER(REST)(__VA_ARGS__)
#define TAIL(a, ...) __VA_ARGS__
#define _TAIL(...) DEFER(TAIL)(__VA_ARGS__)

/* Useful for delays */
#define EMPTY()
/* Useful for generating deliminators */
#define COMMA() ,
/* Useful for consuming macros */
#define VOID(x)

/* Defer macro chains.  Unlikely to need more than 4 */
#define DEFER(m) m EMPTY()
#define DEFER2(m) m EMPTY EMPTY()()
#define DEFER3(m) m EMPTY EMPTY EMPTY()()()
#define DEFER4(m) m EMPTY EMPTY EMPTY EMPTY()()()()

/* Does exactly what you think it does. */
#define EXPAND(...) __VA_ARGS__

/* Concatinator */
#define CAT(a, ...) a ## __VA_ARGS__
#define CONCAT(a, b) _CONCAT(a,b)
#define _CONCAT(a, b) a ## b

/* Hygenic Checker */
#define GENSYM(name) \
  CONCAT(CONCAT(CONCAT(_anon_variable_, name), __LINE__), __COUNTER__)

/* Bool conditional magic below */
#define IS_PROBE(...) SECOND(__VA_ARGS__, 0)
#define PROBE() ~, 1

#define NOT(x) IS_PROBE(CAT(_NOT_, x))
#define _NOT_0 PROBE()

#define BOOL(x) NOT(NOT(x))

#define IF(c) _IF(BOOL(c))
#define _IF(c) CAT(_IF_, c)
#define _IF_0(...)
#define _IF_1(...) __VA_ARGS__

#define IF_ELSE(condition) _IF_ELSE(BOOL(condition))
#define _IF_ELSE(condition) CAT(_IF_ELSE_, condition)

#define _IF_ELSE_1(t, f) t
#define _IF_ELSE_0(t, f) f

#define HAS_ARGS(...) BOOL(FIRST(_END_OF_ARGUMENTS_ __VA_ARGS__)(0))
#define _END_OF_ARGUMENTS_(...) BOOL(FIRST(__VA_ARGS__))

/*--------------------------------------------------------------------------------
  Makes use of C !! magic.
  If we have a pair, NOT_PAIR returns 1.
  If we have greater or fewer than that, NOT_PAIR negates something-not-zero into zero.
  --------------------------------------------------------------------------------*/
#define IS_PAIR(_1, _2, _3, ...) NOT(_3)
#define PAIR(...) BOOL(IS_PAIR(__VA_ARGS__, 0, 1, 1))
#define _PAIR(...) DEFER(PAIR)(__VA_ARGS__)

/* Same logic as IS_PAIR but tells us if we have more than two parameters */
#define IS_GT_TWO(_1, _2, _3, _4, ...) _4
#define GT_TWO(...) BOOL(IS_GT_TWO(1, __VA_ARGS__, 0, 0, 1))
#define _GT_TWO(...) DEFER(GT_TWO)(__VA_ARGS__)


#define MAKE_LIST(...) (__VA_ARGS__)
#define UNDO_LIST(...) __VA_ARGS__

/*--------------------------------------------------------------------------------
  Map functions! Now we can actually do something with all the above insanity.
  Works like you would expect.  Map over the list of values and apply the operator
  and insert the separator to create usable lists.
  --------------------------------------------------------------------------------*/
#define MAP(...)																\
	EVAL(MAP_INNER(__VA_ARGS__))
#define MAP_INNER(op, sep, type, cur_val, ...)	\
  op(type, cur_val)                             \
    IF(HAS_ARGS(__VA_ARGS__))                   \
    (sep()                                      \
     DEFER2(_MAP_INNER)()                       \
     (op, sep, type, ##__VA_ARGS__))
#define _MAP_INNER() MAP_INNER

/* Map over lists of lists.  Calls MAP internally. */
#define MAP_WITH_ID(op,sep, ...)																				\
  IF(HAS_ARGS(__VA_ARGS__))(EVAL(MAP_WITH_ID_INNER(op, sep, ##__VA_ARGS__)))
#define MAP_WITH_ID_INNER(op, sep, cur_val, ...)												\
  DEFER(_MAP_INNER)()(op, sep, FIRST(EXPAND cur_val), TAIL(EXPAND cur_val)) \
    IF(HAS_ARGS(__VA_ARGS__))                                           \
    (sep()                                                              \
     DEFER3(_MAP_WITH_ID_INNER)()                                       \
     (op, sep, FIRST(__VA_ARGS__), TAIL(__VA_ARGS__)))
#define _MAP_WITH_ID_INNER() MAP_WITH_ID_INNER

/* Create a variable with a type.  Useful for declaring a new variable. */
#define MAKE_TYPE_VAR(type, name) type name
/* Create a variable without a type.  Useful for generating parameter lists. */
#define MAKE_NOTYPE_VAR(type, name) name
#define VARS(...) (__VA_ARGS__)

#define VARS_IN(...) DEFER(EXPAND)(__VA_ARGS__)
#define VARS_IN_TYPE(...)                         \
	MAP_WITH_ID(MAKE_TYPE_VAR, COMMA, __VA_ARGS__)
#define VARS_IN_NOTYPE(...)                         \
	MAP_WITH_ID(MAKE_NOTYPE_VAR, COMMA, __VA_ARGS__)


#endif // _DEF_UTIL_H
