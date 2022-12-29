Reviewer 1

I have enjoyed reading about the authors’ software, which I believe has the potential to lead to more robust inferences in observational studies. The following are some relatively minor comments.

Sections 2.2, 2.3, and 2.4– GCOMP, IPTW and AIPTW are commonly used estimators in applications. To improve the readability of the manuscript, I would perhaps shorten and combine sections 2.2, 2.3 and 2.4 into just one section describing existing estimators.

_Response_: Thanks for the suggestion. We have shortened and combined these sections.

Section 3– Unless I am missing something, I would mention that, in addition to the specific moment conditions described by the authors, the estimators need to satisfy also additional and untestable small-bias requirements in order to achieve doubly-robust asymptotic linearity. For example, the second-order terms in Appendix B in Benkeser et al (2017) need to be asymptotically negligible.

_Response_: Thanks, we have added references to the exact regularity conditions as they appear in the original manuscripts. We have also added the sentence 

	> The regularity conditions required to achieve asymptotic linearity include assumptions on the convergence rates of the OR, PS, R-OR, and R-PS to their true counterparts."

Page 6, line starting with “The key theory…”– Should the estimators of the OR and PS satisfy equation (5) rather than (4)?

Page 7, line starting with “Further, they showed…”– Should the estimators of the OR, PS, R-OR satisfy equations (8)-(9) instead of (5)-(6)?

_Response_: Thanks for pointing out. There was some Latex/bookdown referencing gore happening with the original submission. We appreciate the reviewers catching this. It has been corrected in the revision.

Because of the crucial role that the reduced univariate regressions play in obtaining doubly-robust inference, I wonder whether it could be useful to have a simple way of checking the model fit for their estimators. For instance, would it make sense to have functions that can plot the residuals and fitted values against the estimated covariate?

_Response_: The package includes a `plot` method that can be used towards this end. We refer to this method in Section 4.1.1.

Reviewer 2

In this software tutorial, the authors detail and implement state-of-the-art nonparametric estimators of counterfactual means, or the average treatment effect (ATE), from observational studies, together with valid doubly robust statistical inferential tools (confidence intervals and hypothesis tests). Standard inference for doubly robust estimators of the ATE (e.g., AIPW or TML estimators) that employ flexible nuisance estimators is not strictly doubly robust in the same sense as the estimators themselves, i.e., consistency of the estimated ATE only requires consistency of one of two nuisance functions, but this robustness does not apply to asymptotic linearity/normality. That said, valid doubly robust inferential procedures have been proposed by van der Laan (2014) and Benkeser et al. (2017), in which additional nuisance models must be fit and satisfy certain estimating equations. However, these procedures are challenging and involved, thus the drtmle software package fills an important gap by implementing these methods. Beyond these approaches, which extend the usual TML estimator (drtmle function), the authors also include extensions to cross-validated TMLE (cvFolds option, i.e., using sample-splitting), and implement an IPW estimator with valid inference when using a flexible estimator for the propensity score (adaptive_iptw function).

The article does a good job of motivating and describing the underlying methodology well for readers familiar with doubly robust estimators of the ATE, but not with the doubly robust inference literature. Overall, the manuscript is also very clearly written, and the package seems to have a nice design. Specific comments are detailed below.

The workhorse functions (drtmle and adaptive_iptw) allow for missingness in treatment (A) and outcome (Y) values (e.g., with DeltaA and DeltaY control parameters within ‘glm_g’). This bivariate missingness setting doesn’t appear (unless I’m missing something) to have been formally considered in either van der Laan (2014) or Benkeser et al. (2017). I am concerned whether such missingness may lead to non-identifiability of the counterfactual mean—or “treatment-specific marginal mean” as referred to by the authors—functionals, without further assumptions. Are the authors implicitly invoking concrete missingness assumptions under which the implemented approaches target the original functional of interest (e.g., something like Delta_A x Delta_Y independent of (A, Y) given W)? Or are the methods simply standard, or used for convenience? It would enhance the paper if this issue were clarified.

_Response_: Thanks for the suggestion. We have added the formal assumptions needed for identification to that section (copied here): 

if we define $\Delta_A$ and $\Delta_Y$ as indicators that $A$ and $Y$ are observed, respectively, then we can still perform valid inference on $\psi_0$ under the assumptions that (i) $Y(a) \perp \Delta_A \mid W$, (ii) $Y(a) \perp A \mid \Delta_A = 1, W$, and
(iii) $Y(a) \perp \Delta_Y \mid \Delta_A = 1, A = a, W$.

Related to point 1., why is missingness in the covariates (W) not considered? I am aware that this would be a very challenging scenario methodologically (and this is a good enough reason to not deal with it), but I am just curious whether there is something fundamentally different about W in this setting which precludes consideration of it being partially missing. Any guidance on what a practitioner should do if W is partially missing?

_Response_: Thanks for the suggestion. We have added a brief paragraph on missingness in $W$.

In section 4.2 on standard error estimation, the discussion is focused on the AIPW and standard TMLE (i.e., not the DR-TMLE implemented in the package). They say “The ideas generalize immediately to the DR-TML estimator, but the formulas are more complex in this latter case”. Presumably, that means that inference within drtmle is based on the empirical variance of the influence function of the DR-TMLE, but authors should elaborate at least a little. If the formulas are too complex to be presented within the manuscript, readers should be referred to the relevant sections of van der Laan (2014) and/or Benkeser et al. (2017).

_Response_: Thanks for the suggestion. We have restructured this section a bit and refer to the specific location of the formula in Benkeser et al (2017).

Several times, the equations (6) and (7) are referred to as if they were estimating equations, but it seems like they are definitions of new nuisance functions (at the “population” level).

a. On page 7, “Further, they showed that if estimators of the OR, PS, and R-OR satisfy equations (5)–(6)”. Equation (5) is empirical, but 6 is not. Should this say equations (8)—(9)?

b. On page 8, “Now, we have argued that a TML estimator may achieve doubly robust limiting behavior if it additionally satisfies equations (6) and (7)”. Should this say equations (9) and (10)?

c. On Page 8, “iteratively modify initial estimates of the OR and PS via the fluctuateQ() and fluctuateG() functions until equations (5)–(7) are approximately solved”. Should this say equations (8)—(10)?

_Response_: All equation referencing has been corrected. Sorry for the confusion.

Other specific minor comments:

Page 3, “relying upon Neymon orthogonality” should say Neyman?
	
_Response_: Correct. Thanks for catching.

Page 4, “Under additional assumptions, these parameters have interpretations…”. Should at least give a reference here to the causal assumptions.

_Response_: We have added a few short sentences on causal assumptions.

Page 7, equation (9), I think the Q^bar term in the numerator and g_n term in the denominator should each have an asterisk.

_Response_: Yes, thanks. Good catch.

Page 9, “r round(EY0, 2)” looks like unrendered markdown. What is the numerical value?

_Response_: This has been corrected.

Page 18, “methods that compute doubly robust confidence intervals”. Can this be clarified to “Wald-style” confidence intervals? It seems from context that these are just based directly on the standard error estimates from the previous section, but would be good to clarify in opening paragraph.

_Response_: Done
