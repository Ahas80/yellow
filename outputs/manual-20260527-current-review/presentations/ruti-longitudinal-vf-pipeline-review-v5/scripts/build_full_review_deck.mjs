#!/usr/bin/env node

import { buildReviewDeck } from "./longcycler_review_deck_engine.mjs";

function parse(argv) {
  const out = {};
  for (let index = 0; index < argv.length; index += 1) {
    const key = argv[index];
    if (!key.startsWith("--")) throw new Error(`Unexpected positional argument: ${key}`);
    const value = argv[index + 1];
    if (!value || value.startsWith("--")) out[key.slice(2)] = true;
    else { out[key.slice(2)] = value; index += 1; }
  }
  return out;
}

const args = parse(process.argv.slice(2));
buildReviewDeck({
  variant: "v5-full",
  projectRoot: args["project-root"],
  registry: args.registry,
  workspace: args.workspace,
  output: args.output,
  devFixture: args["dev-fixture"],
}).then((result) => console.log(JSON.stringify(result, null, 2))).catch((error) => {
  console.error(error.stack || error.message || String(error));
  process.exit(1);
});

