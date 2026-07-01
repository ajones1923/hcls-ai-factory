# HCLS AI Factory — Skills

This folder holds the **skills** that keep the factory's build aligned with its foundation and
strategy: durable, versioned guidance a developer (or an AI assistant) references while building
the next generation of the platform. Two kinds of thing live here:

- **Invocable skills** — a directory with a `SKILL.md` file (YAML frontmatter + instructions)
  that can be loaded and referenced by name.
- **Reference notes** — plain text/markdown that captures strategy or standing context.

## Layout

```
skills/
├── README.md                          # this file
├── AI Factory Core Pillars.txt        # the 16 foundational pillars (strategic reference)
├── build-housekeeping-standards/      # invocable skill
│   └── SKILL.md                       #   best-practices checklist for every build
└── architecture/                      # (placeholder) architecture skills / notes
```

### `AI Factory Core Pillars.txt`
The 16 pillars the platform is built on — compute, storage, networking, orchestration,
databases, embeddings, messaging, inference serving, observability, security, cost, reliability,
ease of deployment, and GitHub/site presentation. This is the strategic "what the factory is made
of" reference; skills operationalize the housekeeping side of these pillars.

### `build-housekeeping-standards/SKILL.md`
The standing checklist for **how** we build: the pre-commit gate, golden paths for adding an
engine or agent, capability registration (the registry drift-guard), governance wiring, the
dependency baseline, neutrality/secret/path/data hygiene, docs-stay-true, and verify-then-commit
discipline. Reference it whenever adding or modifying an engine, agent, capability, dependency,
test, or doc. It maps back to the Core Pillars (esp. Security, Reliability, Ease of Deployment,
and GitHub Structure).

### `architecture/`
Reserved for architecture-focused skills and notes (diagrams, design decisions, extension
patterns). Empty for now — add a `SKILL.md` here when an architecture skill is written.

## Skill convention

Each invocable skill is its **own directory** containing a `SKILL.md`:

```markdown
---
name: my-skill-name          # kebab-case, matches the directory name
description: >-              # one clear paragraph on WHEN to use it (drives discovery)
  What this skill is for and the situations that should trigger it.
---

# Title

The instructions / checklist body.
```

Keep skills **neutral** (no proprietary vendor or alternate-edition branding — the same rule the
repo's pre-commit guard enforces) so they commit cleanly and stay publishable. Reference real
files by path; keep them actionable and scannable.

## Where skills live and how they're discovered

- **In the repo (`skills/`)** — versioned, reviewable, and shareable. This is the canonical copy
  and what ships with the open-source project.
- **In `.claude/skills/`** — how Claude Code auto-discovers invocable skills in a session.
  `.claude/` is gitignored (local only), so hooking a skill up for auto-load never publishes.

To make a repo skill auto-loadable without duplicating it, symlink it (single source of truth):

```bash
mkdir -p .claude/skills
ln -sfn ../../skills/<skill-name> .claude/skills/<skill-name>
```

`build-housekeeping-standards` is already wired this way. Start a fresh Claude Code session for
newly added skills to be picked up.

## Adding a new skill

1. Create `skills/<skill-name>/SKILL.md` with the frontmatter + body above.
2. Verify it's neutral and passes the repo's checks (`./scripts/install-hooks.sh` then a trial
   `git add` + commit runs the guard).
3. Symlink it into `.claude/skills/` (above) if it should be auto-loaded in-session.
4. Add a one-line entry to the Layout section of this README.
