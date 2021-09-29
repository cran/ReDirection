#usethis::use_testthat()
#usethis::use_test("report")


smat_1 <- matrix(c(1,0,0,0,0,1,0,-1,0,1,0,0,1,0,0,0,-1,1,1,0,0,0,0,0,0,-2,0,-1,-1,1,0,0,0,0,1,0,-1,1,0,0,0,0,0,0,-1,0,1,0,-1,0,0,0,0,0,-1,0,1,0,0,-1,0,-1,0),byrow=TRUE,nrow=9,ncol=7)
testthat::test_that("Algorithm is working.....", {
rates_1 <- calculate_reaction_vector(smat_1,TRUE)  
expect_equal(rates_1,1)
})


testthat::test_that("Algorithm is working.....", {
rates_2 <- calculate_reaction_vector(smat_1)  
expect_equal(rates_2,0)
})


testthat::test_that("Algorithm is working.....", {
rates_3<- calculate_reaction_vector(TRUE)  
expect_equal(rates_3,0)
})


testthat::test_that("Algorithm is working.....", {
rates_4<- calculate_reaction_vector()  
expect_equal(rates_4,0)
})




