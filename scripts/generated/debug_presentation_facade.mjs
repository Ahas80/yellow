import { FileBlob, PresentationFile } from "@oai/artifact-tool";

const starter = "/var/folders/fp/wwwk1rbj70l0k92s92kd2z500000gp/T/codex-presentations/manual-20260710-lecturer-methods/tmp/template-starter.pptx";
const presentation = await PresentationFile.importPptx(await FileBlob.load(starter));
const slide = presentation.slides.items[0];

function describe(label, value) {
  const own = Object.getOwnPropertyNames(value ?? {});
  const proto = value == null ? [] : Object.getOwnPropertyNames(Object.getPrototypeOf(value) ?? {});
  console.log(label, { own, proto });
}

describe("presentation", presentation);
describe("slide", slide);
describe("slide.shapes", slide.shapes);
describe("slide.elements", slide.elements);
console.log("shape count", slide.shapes.items.length);
for (let i = 0; i < slide.shapes.items.length; i += 1) {
  const shape = slide.shapes.items[i];
  console.log("shape", i, {
    id: shape.id,
    name: shape.name,
    position: shape.position,
    text: shape.text?.value ?? shape.text?.text ?? shape.text?.toString?.(),
    textOwn: shape.text ? Object.getOwnPropertyNames(shape.text) : [],
    textProto: shape.text ? Object.getOwnPropertyNames(Object.getPrototypeOf(shape.text) ?? {}) : [],
  });
}
for (const key of Object.getOwnPropertyNames(slide)) {
  const value = slide[key];
  if (value && typeof value === "object") describe(`slide.${key}`, value);
}
