test: test-regression test-unit

test-regression:
	make -C tests/regression test

test-unit:
	make -C tests/unit/Emplode test
