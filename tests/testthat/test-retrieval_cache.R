test_that("Cache is correctly initialized", {
  tmp_dir <- tempdir()
  bfc_nh <- NetworkHub::initialize_NetworkHub(nh_cachedir = tmp_dir)

  expect_s4_class(bfc_nh, "BiocFileCache")
  expect_true(BiocFileCache::bfccount(bfc_nh) == 0)

})
