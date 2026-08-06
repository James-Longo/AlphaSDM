# ee_await_export has two budgets. A backlogged scheduler is worth abandoning
# quickly; a task that is doing real work is not. Conflating them meant a large
# model, which simply takes longer to write, lost the optimisation it most needed.
fake_task <- function(states) {
  i <- 0
  list(
    task = list(
      status = function() { i <<- i + 1; list(state = states[min(i, length(states))]) },
      cancel = function() invisible(NULL)
    ),
    asset_id = "projects/p/assets/a"
  )
}

test_that("a task that never starts is abandoned", {
  expect_error(
    ee_await_export(fake_task(rep("READY", 50)), poll_seconds = 0,
                    max_minutes = 60, max_queue_minutes = 1/600),
    "still queued"
  )
})

test_that("a task that starts is not abandoned for having queued", {
  # The queue budget here is a tenth of a second, and the task runs well past it.
  expect_equal(
    ee_await_export(fake_task(c("READY", rep("RUNNING", 40), "COMPLETED")),
                    poll_seconds = 0, max_minutes = 60, max_queue_minutes = 1/600),
    "projects/p/assets/a"
  )
})

test_that("a running task still stops at the overall limit", {
  expect_error(
    ee_await_export(fake_task(rep("RUNNING", 50)), poll_seconds = 0, max_minutes = 0),
    "exceeded max_minutes"
  )
})

test_that("a failed task reports why", {
  expect_error(
    ee_await_export(fake_task(c("READY", "FAILED")), poll_seconds = 0, max_minutes = 60),
    "FAILED"
  )
})

test_that("no queue limit means queueing alone never ends the wait", {
  expect_equal(
    ee_await_export(fake_task(c(rep("READY", 30), "COMPLETED")), poll_seconds = 0,
                    max_minutes = 60),
    "projects/p/assets/a"
  )
})

test_that("a running task is not cut off by a deadline of our own", {
  # Twice a fixed limit turned Earth Engine finishing the work into this package
  # reporting a failure. There is no overall deadline by default.
  expect_equal(
    ee_await_export(fake_task(c("READY", rep("RUNNING", 500), "COMPLETED")), poll_seconds = 0),
    "projects/p/assets/a"
  )
})

test_that("an unattended run can still impose a ceiling", {
  withr::local_envvar(ALPHASDM_MAX_WAIT_MINUTES = "0.001")
  expect_error(
    ee_await_export(fake_task(rep("RUNNING", 500)), poll_seconds = 0),
    "exceeded max_minutes"
  )
})

test_that("a failed task is reported at once, not waited out", {
  expect_error(
    ee_await_export(fake_task(c("RUNNING", "FAILED")), poll_seconds = 0),
    "FAILED"
  )
})
