context("diagrams")

describe("chunk_ode()", {
  it("separates an ode with leading negative into parts", {
    expect_equal(
      chunk_ode("-ka*Ad+F11*input1"),
      list(structure(c("ka*Ad", "F11*input1"),
                     .Names = c("n", "p")))
    )
  })
  it("separates an ode with no leading symbol into parts", {
    expect_equal(
      chunk_ode("ka*Ad-F11*input1"),
      list(structure(c("ka*Ad", "F11*input1"),
                     .Names = c("p", "n")))
    )
  })
})



