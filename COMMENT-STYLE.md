# COMMENT-STYLE.md

The rules for prose in a repo that adopts this file.  They cover code
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

A letter on a rule number marks a rule inserted at that point in document
order.  It does not make the rule a sub-case of the number it follows.

## Sentences

**1. One thought per sentence.**  A sentence may carry as many clauses as that
thought needs, and joining related clauses usually reads better than splitting
them into short, halting sentences.  The fault this rule names is a sentence
carrying two thoughts that the reader has to pull apart.  Split that one.

**2. Semicolons are rare.**  A semicolon joins two independent clauses, so ask
first whether a comma and a conjunction would carry them, and whether the
second clause wants to be its own sentence.  Keep the semicolon when neither of
those reads as well.

**3. Colons introduce lists and explanations.**  A colon before a real list is
correct, and so is one whose second half explains the clause before it.  What
sits before the colon should be a clause rather than a bare noun phrase used as
a label.

**4. Dashes set off a parenthetical or land a closing clause.**  A pair of
dashes around a parenthetical is acceptable, and real parentheses are usually
better.  A single dash may set off a trailing clause where the pause earns its
keep, but a fragment that wanted to be a sentence is not fixed by hanging it
off a dash.

In a code comment, spell a dash as two hyphens rather than as a real em dash.
Documentation files keep their own practice.

**5. Finish the sentence.**  Detail paragraphs are complete sentences with
verbs.  A `///` brief line or a docstring summary may be a noun phrase.

Keep a brief short.  Brevity is what this asks for, not a line count, so a
brief that still runs past the column limit wraps onto a second line.  Rule
10g covers what a wrapped brief needs from the doc generator.

**5a. A Doxygen brief takes no terminating period.**  That holds whether the
brief reads as a noun phrase or as a full sentence, so every brief in a file
ends the same way.  A `///<` trailing brief follows the same rule.  Detail
paragraphs inside the block punctuate normally.

Another language keeps its own convention.  A Python docstring ends its summary
line with a period, because PEP 257 says so and the linters check it.  What
carries across is that a file settles on one ending, not which ending it is.

## Paragraphs

**6. Keep a paragraph to about five lines.**  Past that, split it into two
paragraphs or cut it down.  A wall of text is not thorough, it is unread.

**7. Put two spaces after a sentence-ending period.**  That holds in code
comments, doc blocks, commit messages, and documentation files alike, markdown
included.  An abbreviation such as "e.g." takes one space, since no sentence
ends there.

Do not bulk-convert a file you are not otherwise editing.  Any line you do
modify comes back with two spaces.

**7a. Keep a line to 96 columns.**  That is a ceiling for code and prose alike.
Comment prose is usually wrapped narrower, and a paragraph you edit keeps the
width it already had rather than being reflowed to the ceiling.  A markdown
file that puts each paragraph on one long line stays that way.

Four things are exempt, because wrapping them costs more than it returns.  A
URL wraps into something the reader cannot click.  A run of column-aligned
declarations or table rows loses the alignment that made it readable.  A string
literal wraps only by splitting it, which puts the text further from what the
program prints.  A pragma or lint suppression that has to sit on one line has
no second line to move to.

**7b. Put one space before a trailing comment in a C-family file.**  That is
the gap between code and the `//` or `///<` that follows it on the same line.
Rule 7 governs a different place, so do not let the two spaces after a period
pull this one to two.  A language with its own rule wins, and PEP 8 asks for two
spaces before an inline `#`.

A run of trailing comments on consecutive declarations is the one exception.
It may pad out to a shared column instead, and then every comment in the run
lines up on that column.  A lone trailing comment is not a run, and neither is
one that follows a statement.

## Content

**8. Do not restate the declaration.**  The reader can see the types, the
parameter names, and the return type.  A comment earns its place by saying
something the code does not.

**9. Do not explain what was rejected.**  Rationale for a road not taken
belongs in a plan document or a ledger, never in a comment or a commit message.

**10. Do not state the obvious.**  "Increment the counter" above a line that
increments the counter is noise.

Cutting has a floor.  Rules 8, 9, and 10 delete what the reader can already
see, never what the reader could not reconstruct.  A comment trimmed until it
is merely shorter than the fact it recorded has been broken, not tightened.

Rule 9 is the one most often overread.  It sends a design alternative that was
weighed and dropped to the ledger.  It does not send away a hazard the reader
is likely to walk into, and a warning against a change that would silently
break something stays next to the code it protects.

**10a. Write an argument once.**  A second site that needs it states the
consequence in a sentence and points at where the argument lives.  Copying the
reasoning means every copy is a place to update when the answer changes.  Point
only at something that travels with the file, so that a file copied into
another program does not end up referring to a document it left behind.

A file written to be copied out must stand alone.  It points at no document
that stays behind, however permanent that document looks from here.  It may
point at a companion file only when it carries that companion along, the way a
header carries the header it includes.

Where the argument lives somewhere that does not travel, the consequence
sentence carries no pointer and stands by itself.  Such a sentence is the only
copy that survives the move, so never thin one into a pointer, and never let a
gotcha live only in a file the copy leaves behind.

One case earns a second copy.  An argument written to head off a suggestion has
to sit where the reader forms the suggestion, because a reader who never follows
the pointer is exactly the one it was written for.  Such a copy carries the
argument alone.  Do not append a line saying the real rationale lives
elsewhere.

## Doxygen

This section assumes Doxygen, and rules 5, 5a, and 7b spell a brief its way.  A
repo on another doc-comment system keeps every rule and reads `///` and `///<`
as whatever that system writes them.  Where that system has its own hard
convention, it wins, and rules 5a and 7b record the two cases that come up.

The Doxygen behavior described below was checked against version 1.18.0.

**10b. The brief goes on `///` lines above the block.**  Write the brief as
`///` immediately above the `/** ... */` block, and let the block hold only the
detail.  Do not use `\brief` inside the block.  A declaration whose brief says
everything needs no block at all.

**10c. Refer to a parameter with `\a`, not `\p`.**  Both mark a word as a
parameter in running text.  `\a` renders it italic and `\p` renders it
monospace.  The italic sets a parameter apart, where monospace blends it into
the identifiers the prose already writes that way.

Neither one checks that the word names a real parameter, so a typo in either is
silent.  Only `\param` checks.  `\p` and the `\c` of rule 10f are one
directive under two names, so this rule bans naming a parameter with it, not
the monospace itself.

**10d. Line the block's stars up with its `/**`.**  A continuation line inside
a `/** ... */` block starts its `*` in the same column as the opening `/**`,
not one column to the right.  The closing `*/` sits in that column too.  A
block indented inside a class indents its stars to match its own `/**`.

**10e. Document an exception with `\exception`, never `\throw` or `\throws`.**
Doxygen renders all three identically.  Committing to one spelling is what
keeps a search for the tag exhaustive.  `\exception` is the noun, which
matches the `\param` and `\return` it sits beside.

**10f. Monospace one token with `\c`, and a phrase with `<code>`.**  `\c` marks
only the word that follows it, so `\c operator new` sets `operator` in
monospace and leaves `new` in the running text.  Doxygen reports nothing when
that happens, the same way it says nothing about a misspelled `\a`.

Reach for `<code>...</code>` when the span runs to more than one token, and for
nothing else.  Two spellings with a rule between them keep a search for either
exhaustive.

**10g. Set `MULTILINE_CPP_IS_BRIEF = YES`.**  Under the default of `NO`, only a
one-line `///` block counts as a brief.  A brief that wraps onto a second line
then produces no brief at all, and its text falls through into the detailed
description, where it runs together with the first detail paragraph.  Doxygen
reports nothing.

A repo that generates no documentation still sets it, in whatever scratch
configuration is used to check a file.  The wrapped briefs rule 5 allows are
correct only under this setting.

**11. Use `\copydoc` when the whole block transfers.**  It copies the `\param`
and `\exception` lines along with the text, so it only works when the parameter
names match.  An overload that names them differently draws an "is not found in
the argument list" warning, and a block with its own `\exception` lines ends up
carrying each one twice.

**12. Use `\copybrief` plus one specific paragraph when they do not.**  Two
overloads that spell the same parameter differently are the usual case.  One
taking `data` and `len` cannot `\copydoc` from one taking `src`, because the
imported `\param` names would not match the signature.

**13. Verify that a reference resolves.**  Doxygen cannot tell two overloads
apart by a `requires` clause alone.  A reference that silently resolves to the
wrong one is worse than no reference.

## Wording

**14. American spellings.**  Watch `-our`, `-ise`, `-re`, and "judgment".
Never "programme", in any sense.

**15. Use the serial comma** before the final conjunction.

**16. Never write "load-bearing".**  Name the actual dependency instead.
