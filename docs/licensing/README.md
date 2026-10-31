# License proposal — pending rights-holder approval

**Status: proposal only. This directory does not grant a license to the existing MATLAB code or change the repository's current licensing status.** The candidate license is quoted as review material, and the README patch has not been applied.

## Proposed choice

Use the unmodified BSD 3-Clause License for the original MATLAB source files and README, subject to confirmation of ownership and authorization. MIT is also a reasonable permissive alternative. Choose GPL-3.0 instead only if the rights holders want covered derivatives distributed to others to retain GPL source-availability and licensing obligations.

| Topic | MIT | BSD-3-Clause | GPL-3.0 |
| --- | --- | --- | --- |
| Academic paper citation | No mandatory paper-citation term in the standard license | Same | Same |
| Modify the implementation | Allowed | Allowed | Allowed; private changes generally need not be published |
| Source redistribution | Retain copyright and permission notices | Retain copyright, conditions and disclaimer | Preserve required notices and license text; covered derivative works must satisfy GPL terms |
| Binary redistribution | Preserve applicable notices; no source-disclosure requirement | Reproduce required notices in accompanying materials; no source-disclosure requirement | Provide Corresponding Source or a compliant means of obtaining it under section 6 |
| Commercial use | Allowed, including compliant proprietary integration | Same | Allowed, but a covered combined derivative cannot be distributed solely under proprietary terms |

BSD-3-Clause additionally expressly restricts using copyright holders' or contributors' names to endorse or promote derived products without specific prior written permission. This is distinct from ordinary academic citation.

## Review materials

- [Candidate BSD-3-Clause text](LICENSE.bsd-3-clause.proposed.md): a quoted template, with the copyright holders deliberately unfilled.
- [Proposed README patch](README.proposed.patch): preserves the existing BibTeX and experiment instructions; separates citation/reproducibility requests from license conditions and clarifies the MATLAB/Qhull boundary.

The patch was prepared against commit `fe12e2734899716fff50861e9ebf1273a98a25e0`; the README blob is `0ead30088511a3fd644519f85b75fb795826e24c`.

## Dependency findings

The reviewed tree contains 32 MATLAB files and no bundled third-party package, MEX binary, or submodule. MATLAB functions such as `pcg`, `eigs`, and `delaunayTriangulation` are used. MathWorks documents that `convhulln` and `delaunayn` use Qhull; this repository does not distribute Qhull source or binaries. The local `sym` function at the end of `mg_setup.m` is numerical symmetrization, not evidence of a Symbolic Math Toolbox dependency.

These are static inspection findings, not proof of originality or a runtime dependency certification. No MATLAB experiments were run for this documentation change. Licensing this repository does not relicense MATLAB, its supplied software, or the linked preprint. GPL use with a proprietary interpreter and any library bindings must be assessed according to the actual combination rather than assumed universally prohibited or universally compatible.

## Decisions required before publishing an operative license

- [ ] Identify the actual copyright holders for the code and README; paper authorship and repository administration do not by themselves establish ownership.
- [ ] Confirm permission from all necessary contributors, employers/institutions, or other relevant rights holders.
- [ ] Check copied/adapted snippets and any third-party notices; absence of a notice is not proof of original authorship.
- [ ] Confirm acceptance of BSD-3-Clause, including commercial use and compliant proprietary integration without mandatory upstream contribution.
- [ ] Confirm copyright names, years, scope, and any relevant patent or contractual arrangements.

Once these decisions are recorded, replace the placeholder in the candidate text, place the standard text in a root `LICENSE`, and apply the README patch in a separate explicitly authorized change. Do not publish the candidate placeholder as an operative license. The MIT license of a separate academic-homepage repository does not authorize this code.

## Primary references

- [MIT license — Open Source Initiative](https://opensource.org/license/mit)
- [BSD-3-Clause — Open Source Initiative](https://opensource.org/license/bsd-3-clause)
- [GPL version 3 — Open Source Initiative](https://opensource.org/license/gpl-3.0)
- [GNU license FAQ — Free Software Foundation](https://www.gnu.org/licenses/gpl-faq.en.html#InterpreterIncompat)
- [Licensing a repository — GitHub Docs](https://docs.github.com/en/repositories/managing-your-repositorys-settings-and-features/customizing-your-repository/licensing-a-repository)
- [convhulln — MathWorks](https://www.mathworks.com/help/matlab/ref/convhulln.html)
- [delaunayn — MathWorks](https://www.mathworks.com/help/matlab/ref/delaunayn.html)
