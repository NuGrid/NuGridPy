# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line.
SPHINXOPTS    =
SPHINXBUILD   = sphinx-build
SPHINXPROJ    = PPM
SOURCEDIR     = source
BUILDDIR      = build

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: help Makefile

gh-pages:
	$(shell bash -c 'read -s -p "THIS PROCESS PUSHES YOUR LOCAL REPO TO THE MASTER ON GITHUB. Would you like to continue? [enter] = continue, [^C] = exit: " pwd; echo $$pwd')
	make clean
	make html
	git add -A && git commit -m "building and pushing docs" && git push origin master
	cd .. && git checkout gh-pages && touch .nojekyll && rm -rf _sources _static _images _autosummary && git checkout master docs/build/html && mv ./docs/build/html/* ./ && rm -rf ./docs && git add -A && git commit -m "Generated gh-pages for `git log master -1 --pretty=short --abbrev-commit`" && git push origin gh-pages ; git checkout master
    
# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%: Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)