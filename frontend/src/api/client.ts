export type ExpressionMatrix = Record<string, Record<string, number>>
export type GeneSets = Record<string, string[]>

// API base URL from env (Vite exposes only VITE_ prefixed vars). Empty string means
// use relative URLs via Vite proxy (/api -> :8000). Do not hardcode localhost.
const API_BASE: string = (import.meta.env.VITE_API_URL as string | undefined)?.replace(/\/$/, '') ?? ''

export function getApiBase(): string {
  return API_BASE
}

export async function scoreBoth(
  expression: ExpressionMatrix,
  geneSets: GeneSets,
  alpha = 0.25,
  baseUrl: string = API_BASE
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
  baseUrl: string = API_BASE
): Promise<any> {
  const res = await fetch(`${baseUrl}/api/v1/differential`, {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify({ scores, groups }),
  })
  if (!res.ok) throw new Error(await res.text())
  return res.json()
}

export async function health(baseUrl: string = API_BASE): Promise<any> {
  const r = await fetch(`${baseUrl}/health`)
  return r.json()
}
