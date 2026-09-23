# ee_await_tasks() waits for every task and reports failures instead of
# stopping, so each caller decides what a failure means (the Drive export
# treats an all-water cell as empty, a table export stops).
fake_task <- function(states, error = NULL) {
  i <- 0
  list(status = function() {
    i <<- i + 1
    st <- states[min(i, length(states))]
    c(list(state = st), if (st %in% c("FAILED", "CANCELLED") && !is.null(error))
      list(error_message = error))
  })
}

test_that("tasks that complete, after queueing, report no failures", {
  tasks <- list(a = fake_task(c("READY", "RUNNING", "COMPLETED")),
                b = fake_task("COMPLETED"))
  out <- suppressMessages(ee_await_tasks(tasks, poll_seconds = 0))
  expect_length(out, 0)
})

test_that("a failed task is reported with its message and the rest still finish", {
  tasks <- list(ok  = fake_task(c("RUNNING", "COMPLETED")),
                bad = fake_task(c("RUNNING", "FAILED"), error = "User memory limit exceeded."))
  out <- suppressMessages(ee_await_tasks(tasks, poll_seconds = 0))
  expect_identical(names(out), "bad")
  expect_match(out[["bad"]], "memory limit")
})

test_that("a cancelled task counts as a failure", {
  out <- suppressMessages(ee_await_tasks(list(t = fake_task("CANCEL_REQUESTED")),
                                         poll_seconds = 0))
  expect_identical(unname(out), "CANCELLED")
})
