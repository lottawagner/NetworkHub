test_that("Cache is correctly initialized", {
  tmp_dir <- tempdir()
  bfc_nh <- NetworkHub::initialize_NetworkHub(nh_cachedir = tmp_dir)

  expect_s4_class(bfc_nh, "BiocFileCache")
  expect_true(BiocFileCache::bfccount(bfc_nh) == 0)

  bfc_nh

  cpdb_url <- urlmaker_cpdb(species = "human")
  test_rpath <- BiocFileCache::bfcadd(x = bfc_nh, rname = "cpdbtest", cpdb_url)

  expect_true(BiocFileCache::bfccount(bfc_nh) == 1)
  expect_type(test_rpath, "character")
})
