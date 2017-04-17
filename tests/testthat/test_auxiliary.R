##########################################
# tests for auxillary functions in...
# 1. auxiliaries_io.R
# 2. auxiliaries_string.R
# 3. auxiliaries_vector.R
##########################################

context("auxiliaries_io")

test_that("'fileparts' parses appropriately", {
  fp_out <- fileparts("a/b/c/d.R")
  expect_equal(fp_out$pathname, "a/b/c")
  expect_equal(fp_out$filename, "d")
  expect_equal(fp_out$fileext, ".R")

  fp_out <- fileparts("name")
  expect_equal(fp_out$pathname, ".")
  expect_equal(fp_out$filename, "name")
  expect_equal(fp_out$fileext, "")

  fp_out <- fileparts("name.R")
  expect_equal(fp_out$pathname, ".")
  expect_equal(fp_out$filename, "name")
  expect_equal(fp_out$fileext, ".R")

  fp_out <- fileparts(".R")
  expect_equal(fp_out$pathname, ".")
  expect_equal(fp_out$filename, "")
  expect_equal(fp_out$fileext, ".R")
})

# testing on most input-output functions has been omitted for now.

context("auxiliaries_string.R")

test_that("'strtrimM' removes leading/trailing whitespace", {
  cond1 <- c("a b c d e")
  cond2 <- c("abcde")
  expect_equal(strtrimM("a b c d e"), cond1)
  expect_equal(strtrimM(intToUtf8(c(9, 10, 11, 12, 13, 32, 97,  32,  98,  32,  99,  32, 100,  32, 101, 32, 13))),
               cond1)
  expect_equal(strtrimM(" a b c d e    "), cond1)
  expect_equal(strtrimM("abcde"), cond2)
  expect_equal(strtrimM("abcde   "), cond2)
  expect_equal(strtrimM("    abcde"), cond2)
  expect_equal(strtrimM("   abcde   "), cond2)
  expect_equal(strtrimM("   ab c de   "), c("ab c de"))
})

test_that("'strrepM' replaces elements within a string appropriatley", {
  expect_equal(strrepM("a b c d e"," ",""), c("abcde"))
  expect_equal(strrepM("Hello World!","World","Kitty"), "Hello Kitty!")
})

test_that("'strremWhite' removes all whitespace in a string", {
  cond <- c("abcde")
  expect_equal(strremWhite("abcde"), cond)
  expect_equal(strremWhite(" abcde  "), cond)
  expect_equal(strremWhite("a b c d e"), cond)
  expect_equal(strremWhite(" a b c d e    "), cond)
  expect_equal(strremWhite(intToUtf8(c(9, 10, 11, 12, 13, 32, 97,  32,  98,  32,  99,  32, 100,  32, 101, 32, 13))),
               cond)
})

test_that("'strmatch' returns correct indicies", {
  expect_equal(strmatch(c("hello","test","a","e","hello"),"a"), 3)
  expect_equal(strmatch(c("hello","test","a","e","hello"),"hello"), c(1,5))
  expect_equal(strmatch(c("hello","test","a","e","hello"),"ello"), NULL)
  expect_equal(strmatch(c("hello","test","a","e","hello"),"Hello"), NULL)
  expect_equal(strmatch(c("hello","test","a","e","hello"),"xyz"), NULL)
  expect_equal(strmatch(c(1,2,3,4,8,3),3),c(3,6))
})

test_that("'strexplode' splits string according to a separator", {
  expect_equal(strexplode(c("Hello,,Test,,X,,Z")), c("Hello", "", "Test", "", "X", "", "Z"))
  expect_equal(strexplode(c("Hello,,Test,,X,,Z"),","), c("Hello", "", "Test", "", "X", "", "Z"))
  expect_equal(strexplode(c("Hello,,Test,,X,,Z"),"\\|"), c("Hello,,Test,,X,,Z"))
  expect_equal(strexplode(c("Hello,,Test"),"|"),
               c("H", "e", "l", "l", "o", ",", ",", "T", "e", "s", "t"))  # intended behavior?
  expect_equal(strexplode(c("Hello,,Test,,X,,Z"),",,"), c("Hello", "Test", "X", "Z"))
  expect_equal(strexplode(c("Hello,,Test"),","), c("Hello", "", "Test"))  # intended behaviour?
})

test_that("'strexplodePC' treats parentheses appropriately", {
  expect_equal(strexplodePC("(a,b),(c,d)"), c("(a,b)","(c,d)"))
  expect_equal(strexplodePC("[a,b],[c,d]"), c("[a", "b]", "[c", "d]"))
  expect_equal(strexplodePC("{a,b},{c,d}"), c("{a", "b}", "{c", "d}"))
  expect_equal(strexplodePC("((a,b),(c,d))"), c("((a,b),(c,d))"))
  expect_equal(strexplodePC("(a,b),(c,d))"), c("(a,b)", "(c,d))"))
  expect_equal(strexplodePC("((a,b),(c,d)"), c("((a,b),(c,d)"))  # intended behavior?
  expect_equal(strexplodePC("(a;b);(c;d)",";"), c("(a;b)", "(c;d)"))
  expect_equal(strexplodePC("(a;b),(c;d)",";"), c("(a;b),(c;d)"))
  expect_equal(strexplodePC("(a,b)|(c,d)","|"), c("(a,b)", "(c,d)"))
})

test_that("'strlocateall' returns approprate list of indices", {
  cond1 <- strlocateall("abcdefgabcdefg","a")
  cond2 <- strlocateall("abcdefgabcdefg","ab")
  cond3 <- strlocateall("abcdefgabcdefg","abcdefgabcdefg")

  expect_equal(cond1$start, c(1,8))
  expect_equal(cond1$end, c(1,8))

  expect_equal(cond2$start, c(1, 8))
  expect_equal(cond2$end, c(2, 9))

  expect_equal(cond3$start, 1)
  expect_equal(cond3$end, 14)

  expect_null(strlocateall("abcdefgabcdefg","x")$start)
  expect_null(strlocateall("abcdefgabcdefg","x")$end)
})

context("auxiliaries_vector.R")

test_that("is_numeric returns a numeric vector (with NA if necessary)", {
  expect_equal(as_numeric(c(1,2,3)), c(1,2,3))
  expect_equal(as_numeric(c('1',2,3)), c(1,2,3))
  expect_equal(as_numeric(c('1','a',3)), c(1, NA, 3))
  expect_equal(as_numeric(c('1','2',NA)), c(1, 2, NA))
  expect_equal(as_numeric(factor(c(1, "4", 3))), c(1, 4, 3))
  expect_equal(as_numeric(factor(c(1, "4", NA))), c(1, 4, NA))
  expect_equal(as_numeric(c('1',',',3)), c(1, NA, 3))
})

test_that("'isnumericVector' detects numeric elements appropriately", {
  expect_true(isnumericVector(c(1,2,3)))
  expect_true(isnumericVector(c('1',2,3)))
  expect_true(isnumericVector(factor(c(1, "4", 3))))

  expect_false(isnumericVector(c('1','a',3)))
  expect_false(isnumericVector(c('1','2',NA)))
})
