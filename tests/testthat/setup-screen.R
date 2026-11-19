# The app now opens the new main screen by default (.pt_screen(), R/app-screen-ui.R). The
# older tests build and drive the classic composer through .pt_composer_app() (its
# inputs, its state, the names in environment(app$serverFuncSource())), so a test run
# selects the classic screen, here, until the run ends. A setup file, not a helper:
# pkgload::load_all() sources helpers too, which would switch a developer's own session
# to the classic screen. A test of the new screen opts in locally:
# withr::local_options(PlotTrellis.screen = "new").
withr::local_options(list(PlotTrellis.screen = "classic"), .local_envir = testthat::teardown_env())
