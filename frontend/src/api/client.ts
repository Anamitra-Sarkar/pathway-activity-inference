export type ExpressionMatrix = Record<string, Record<string, number>>
export type GeneSets = Record<string, string[]>

export async function scoreBoth(
  expression: ExpressionMatrix,
  geneSets: GeneSets,
  alpha = 0.25,
  baseUrl = ''
): Promise<any> {
  const res = await fetch(`${baseUrl}/api/v1/score/both`, {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify({ expression, gene_sets: geneSets, alpha }),
  })
  if (!res.ok) throw new Error(await res.text())
  return res.json()
}

export async function differential(
  scores: Record<string, Record<string, number>>,
  groups: Record<string, string>,
  baseUrl = ''
): Promise<any> {
  const res = await fetch(`${baseUrl}/api/v1/differential`, {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify({ scores, groups }),
  })
  if (!res.ok) throw new Error(await res.text())
  return res.json()
}

export async function health(baseUrl = ''): Promise<any> {
  const r = await fetch(`${baseUrl}/health`)
  return r.json()
}
