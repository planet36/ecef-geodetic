# COMMENT-STYLE.md

This file sets the rules for prose in a repo that adopts it.  They cover code
comments, doc blocks, and commit messages.  Documentation files follow them
too, though a table or a reference list may need its own shape.

These rules exist because one habit produced almost every prose complaint
raised against the code they were written for.  That habit is packing several
ideas into one sentence and welding them together with punctuation instead of
writing several sentences.  The cure is not a ban on complex sentences.
Clauses that belong to one thought are better joined than chopped apart, and
the rules below are about where that line falls.

Write the rules into the first draft.  A comment that needs a rewrite pass to
satisfy them was written wrong.

A letter on a rule number marks a sub-case of the rule with that number.

## Sentences

**1. One thought per sentence.**  A sentence may carry as many clauses as that
thought needs, and joining related clauses usually reads better than splitting
them into short, halting sentences.  The fault this rule names is a sentence
carrying two thoughts that the reader has to pull apart.  Split that one.

**2. Semicolons are rare.**  A semicolon joins two independent clauses, so ask
first whether a comma and a conjunction would carry them, and whether the
second clause wants to be its own sentence.  Keep the semicolon when neither of
those reads as well.

**3. Colons introduce rather than join.**  A colon before a real list is
correct, and so is one before a block set off below it, such as a code sample
or a quotation.  So is one after a conventional label, such as "Usage:".  Do not
use a colon to join a clause to its explanation, as in "X must match the word
size: the largest power of 2 bits".  Write "which is", a conjunction, or a
second sentence instead, and do not move the join to a dash.  Apart from such a
label, what sits before the colon should be a clause rather than a bare noun
phrase.  A bold run-in label ends with a period rather than a colon, as each
rule's title here does, and the text after it is a full sentence.

**4. Dashes set off a parenthetical or land a closing clause.**  A pair of
dashes around a parenthetical is acceptable, and real parentheses are usually
better.  A single dash may set off a trailing clause where the pause earns its
keep, but a fragment that wanted to be a sentence is not fixed by hanging it
off a dash.

In a code comment or a commit message, spell a dash as two hyphens rather than
as a real em dash.  Documentation files keep their own practice.

**5. Finish the sentence.**  Detail paragraphs are complete sentences with
verbs.  A summary line is the exception.  Rule 21 makes a commit subject an
imperative phrase, rule 22 lets a Doxygen brief be a noun phrase, and a
docstring summary follows its language's convention.

## Paragraphs

**6. Keep a paragraph to one topic.**  A wall of text is not thorough — it is
unread.  Settle the sentences under rule 1 first, because splitting a paragraph
of run-on sentences only spreads them out.

Look closer at a paragraph past about 100 words or 5 sentences, counting
rendered words rather than source lines, since wrapping varies from file to
file.  Neither number is an automatic split.  Split where the topic, the step,
or the audience changes, and open each new paragraph with the sentence that
states its point.

A sequence of steps or parallel items written as prose is a list, so make it
one.  A list item that is a full sentence ends with a period, and a fragment
does not.  One list does not mix the two.

Keep an example or a caveat with the claim it qualifies.  A single idea stays
in one paragraph even when it runs long, and if it runs too long, cut it rather
than split it.  A split changes paragraph breaks, not wording.

## Punctuation

The marks that join clauses into one sentence are under Sentences, and this
section holds the rest.

**7. Put two spaces after a sentence-ending period.**  That holds in code
comments, doc blocks, commit messages, and documentation files alike, Markdown
included.  An abbreviation such as "vs." takes one space, since no sentence
ends there.

Do not bulk-convert a file you are not otherwise editing.  Any line you do
modify comes back with two spaces.

**8. Put punctuation inside quotation marks only when it belongs to the quoted
text.**  American spelling does not bring with it the American habit of
putting every period and comma inside.  Quoted text in code and commit messages
is usually literal, such as a commit subject, an error message, or a command,
and a mark added inside the quotes reads as part of it.

**9. Use the serial comma.**  Put a comma before the final conjunction in a
list of three or more.

**10. Follow "e.g." and "i.e." with a comma.**  Each one stands in for "for
example" or "that is", which would take a comma in the same place.

## Layout

**11. Keep a line to 96 columns.**  That is a ceiling for code and prose alike.
Comment prose is usually wrapped narrower, and a paragraph you edit keeps the
width it already had rather than being reflowed to the ceiling.  A Markdown
file that puts each paragraph on one long line stays that way.

Four more things are exempt, because wrapping them costs more than it returns.
A URL wraps into something the reader cannot click.  A run of column-aligned
declarations or table rows loses the alignment that made it readable.  A string
literal wraps only by splitting it, which puts the text further from what the
program prints.  A pragma or lint suppression that has to sit on one line has
no second line to move to.

**12. Put one space before a trailing comment in a C-family file.**  That is
the gap between code and the `//` or `///<` that follows it on the same line.
Rule 7 governs a different place, so do not let the two spaces after a period
pull this one to two.  A language with its own rule wins, and PEP 8 asks for at
least two spaces before an inline `#`.

A run of trailing comments on consecutive declarations is the one exception.
It may pad out to a shared column instead, and then every comment in the run
lines up on that column.  A lone trailing comment is not a run, and neither is
one that follows a statement.

## Content

**13. Do not restate the declaration.**  The reader can see the types, the
parameter names, and the return type.  A comment earns its place by saying
something the code does not.

**14. Do not explain what was rejected.**  Rationale for a road not taken
belongs in a plan document or a decision ledger, never in a comment or a
commit message.

**15. Do not state the obvious.**  "Increment the counter" above a line that
increments the counter is noise.

Cutting has a floor.  Rules 13 and 15 delete what the reader can already see,
never what the reader could not reconstruct.  A comment trimmed until it is
merely shorter than the fact it recorded has been broken, not tightened.

Rule 14 is the one most often overread.  It sends a design alternative that
was weighed and dropped to the ledger.  It does not send away a hazard the
reader is likely to walk into, and a warning against a change that would
silently break something stays next to the code it protects.  Nor does it send
away an argument that heads off a likely suggestion, which rule 16 keeps where
the suggestion forms.

**16. Write an argument once.**  A second site that needs it states the
consequence in a sentence and points at where the argument lives.  Copying the
reasoning means every copy is a place to update when the answer changes.  Point
only at something that goes along when the file is copied into another program.

A document that stays behind is off limits, however permanent it looks from
here.  A file may point at a companion file only when it carries that companion
along, the way a header carries the header it includes.

Where the argument lives somewhere that does not travel, the consequence
sentence carries no pointer and stands by itself.  Such a sentence is the only
copy that survives the move, so never thin one into a pointer, and never let a
gotcha live only in a file the copy leaves behind.

One case earns a second copy.  An argument written to head off a suggestion has
to sit where the reader forms the suggestion, because a reader who never follows
the pointer is exactly the one it was written for.  Such a copy carries the
argument alone.  Do not append a line saying the real rationale lives
elsewhere.

## Wording

**17. Use American spellings.**  Watch for the British endings `-our`, `-ise`,
`-yse`, and `-re`, and for a doubled `l` before a suffix, as in "travelled".
These spellings are forbidden in any sense: "programme", "judgement",
"licence", "defence", and "grey".  A spelling that is code, such as an
identifier, stays as the code spells it.

**18. Never write "load-bearing".**  Name the actual dependency instead.

**19. Write numbers in one format.**  A number of four digits stands plain, and
one of five or more digits takes a comma before each group of three, as in 6200
and 14,400.  A quantity the reader might check or tune reads best as numerals,
as in 2 iterations, 96 columns, and 100 words, and a small everyday count reads
best as a word, as in two spaces.  A unit symbol takes a space after the number,
as in 10 nm, but a degree or percent sign takes none, as in 90° and 5%.  A
number that is code, such as a command-line argument or a literal, stays as the
code spells it.

**20. Write a heading in sentence case.**  Capitalize its first word and
lowercase the rest, as in "Benchmark results".  A name, an acronym, or a word
that is code keeps its own spelling wherever it falls, so "glibc" stays
lowercase even as the first word.  A heading that quotes a title, such as a
section of a standard, keeps the title's own case.

## Commit messages

**21. A commit subject is an imperative phrase with no terminating period.**
It names what the commit does, as in "Add a README".  It is not past tense, as
in "Added a README", and not a description, as in "Adds a README".  It is the
commit's brief, so it ends without a period, as a Doxygen brief does.

The body below it is detail and punctuates normally.  It may tell what the
code did before in the past tense, and what it does now in the present.  An
editor that recognizes a commit message wraps its body at 72 columns, so a body
reads best at that width.

## Doxygen

This section assumes Doxygen, and so does the `///<` in rule 12.  A repo on
another doc-comment system keeps every rule and reads `///` and `///<` as
whatever that system writes them.  Where that system has its own hard
convention, it wins, and rules 12 and 22a record the two cases that arise.

The Doxygen behavior described below was checked against version 1.18.0.

**22. Write a short brief on `///` lines above the block.**  Put the brief on
`///` immediately above the `/** ... */` block, and let the block hold only the
detail.  Do not use `\brief` inside the block.  A declaration whose brief says
everything needs no block at all.

A brief may be a noun phrase.  Keep it short.  Brevity is what this asks for,
not a line count, so a brief that still runs past the column limit wraps onto a
second line.  Rule 22b covers what a wrapped brief needs from the doc
generator.

**22a. A Doxygen brief takes no terminating period.**  That holds whether the
brief reads as a noun phrase or as a full sentence, so every brief in a file
ends the same way.  A `///<` trailing brief follows the same rule.  Detail
paragraphs inside the block punctuate normally.

Another language keeps its own convention.  A Python docstring ends its summary
line with a period, because PEP 257 says so.  What carries across is that a
file settles on one ending, not which ending it is.

**22b. Set `MULTILINE_CPP_IS_BRIEF = YES`.**  Under the default of `NO`, only a
one-line `///` block counts as a brief.  A brief that wraps onto a second line
then produces no brief at all, and its text falls through into the detailed
description, where it runs together with the first detail paragraph.  Doxygen
reports nothing.

A repo that generates no documentation still sets it, in whatever scratch
configuration is used to check a file.  The wrapped briefs rule 22 allows are
correct only under this setting.

**23. Use `\a` for a parameter and `\c` for a code token, never `\p`.**  `\a`
and `\p` both mark a word as a parameter in running text.  `\a` renders it
italic and `\p` renders it monospace.  The italic sets a parameter apart, where
monospace blends it into the identifiers the prose already writes that way.

Neither one checks that the word names a real parameter, so a typo in either is
silent.  Only `\param` checks.  `\p` is another name for the `\c` of rule 24,
so for code it would only add a second spelling to search for.

**24. Monospace one token with `\c`, and a phrase with `<code>`.**  `\c` marks
only the word that follows it, so `\c operator new` sets `operator` in
monospace and leaves `new` in the running text.  Doxygen reports nothing when
that happens, the same way it says nothing about a misspelled `\a` parameter.

Reach for `<code>...</code>` when the span runs to more than one token, and for
nothing else.  Two spellings with a rule between them keep a search for either
exhaustive, so do not add Markdown backticks as a third.

**25. Line the block's stars up with its `/**`.**  A continuation line inside
a `/** ... */` block starts its `*` in the same column as the opening `/**`,
not one column to the right.  The closing `*/` sits in that column too.  A
block indented inside a class indents its stars to match its own `/**`.

**26. Document an exception with `\exception`, never `\throw` or `\throws`.**
Doxygen renders all three identically.  Committing to one spelling is what
keeps a search for the tag exhaustive.  `\exception` is the noun, which
matches the `\param` and `\return` it sits beside.

**27. Use `\copydoc` when the whole block transfers.**  It copies the `\param`
and `\exception` lines along with the text, so it only works when the parameter
names match.  An overload that names them differently draws an "is not found in
the argument list" warning, and a block with its own `\exception` lines ends up
carrying each one twice.

Matching names is the only part Doxygen checks.  A `\pre` or a postcondition
that holds for only some of the overloads is copied into the rest, and so is an
`\a` naming a parameter that some of them lack.  Doxygen reports neither.  A
target block therefore holds only what is true of every overload that copies
it.  Where the contracts differ, the whole block does not transfer, and rule
27a applies.

**27a. Use `\copybrief` plus one specific paragraph when it does not.**  Two
overloads that spell the same parameter differently are the usual case.  One
taking `data` and `len` cannot `\copydoc` from one taking `src`, because the
imported `\param` names would not match the signature.

**28. Verify that a reference resolves.**  Doxygen cannot tell two overloads
apart by a `requires` clause alone.  A reference that silently resolves to the
wrong one is worse than no reference.
