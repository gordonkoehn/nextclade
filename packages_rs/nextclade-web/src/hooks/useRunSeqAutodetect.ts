import { useRouter } from 'next/router'
import { useRecoilCallback } from 'recoil'
import { ErrorInternal } from 'src/helpers/ErrorInternal'
import { axiosFetch } from 'src/io/axiosFetch'
import { autodetectResultAtom, autodetectResultsAtom, minimizerIndexAtom } from 'src/state/autodetect.state'
import { minimizerIndexVersionAtom } from 'src/state/dataset.state'
import { MinimizerIndexJson, MinimizerSearchRecord } from 'src/types'
import { qrySeqInputsStorageAtom } from 'src/state/inputs.state'
import { getQueryFasta } from 'src/workers/launchAnalysis'
import type { Subscription } from 'observable-fns'
import { NextcladeSeqAutodetectWasmWorker } from 'src/workers/nextcladeAutodetect.worker'
import { spawn } from 'src/workers/spawn'

export function useRunSeqAutodetect() {
  const router = useRouter()

  return useRecoilCallback(
    ({ set, reset, snapshot: { getPromise } }) =>
      () => {
        reset(minimizerIndexAtom)
        reset(autodetectResultsAtom)

        void router.push('/autodetect', '/autodetect') // eslint-disable-line no-void

        function onResult(res: MinimizerSearchRecord) {
          set(autodetectResultAtom(res.fastaRecord.index), res)
        }

        Promise.all([getPromise(qrySeqInputsStorageAtom), getPromise(minimizerIndexVersionAtom)])
          .then(async ([qrySeqInputs, minimizerIndexVersion]) => {
            if (!minimizerIndexVersion) {
              throw new ErrorInternal('Tried to run minimizer search without minimizer index available')
            }
            const fasta = await getQueryFasta(qrySeqInputs)
            const minimizerIndex: MinimizerIndexJson = await axiosFetch(minimizerIndexVersion.path)
            set(minimizerIndexAtom, minimizerIndex)
            return runAutodetect(fasta, minimizerIndex, onResult)
          })
          .catch((error) => {
            throw error
          })
      },
    [router],
  )
}

async function runAutodetect(
  fasta: string,
  minimizerIndex: MinimizerIndexJson,
  onResult: (res: MinimizerSearchRecord) => void,
) {
  const worker = await SeqAutodetectWasmWorker.create(minimizerIndex)
  await worker.autodetect(fasta, { onResult })
  await worker.destroy()
}

export class SeqAutodetectWasmWorker {
  private thread!: NextcladeSeqAutodetectWasmWorker
  private subscription?: Subscription<MinimizerSearchRecord>

  private constructor() {}

  static async create(minimizerIndex: MinimizerIndexJson) {
    const self = new SeqAutodetectWasmWorker()
    await self.init(minimizerIndex)
    return self
  }

  async init(minimizerIndex: MinimizerIndexJson) {
    this.thread = await spawn<NextcladeSeqAutodetectWasmWorker>(
      new Worker(new URL('src/workers/nextcladeAutodetect.worker.ts', import.meta.url), {
        name: 'nextcladeAutodetectWorker',
      }),
    )

    await this.thread.create(minimizerIndex)
  }

  async autodetect(
    fastaStr: string,
    {
      onResult,
      onError,
      onComplete,
    }: {
      onResult: (r: MinimizerSearchRecord) => void
      onError?: (error: Error) => void
      onComplete?: () => void
    },
  ) {
    this.subscription = this.thread.values().subscribe(onResult, onError, onComplete)
    await this.thread.autodetect(fastaStr)
  }

  async destroy() {
    await this.subscription?.unsubscribe() // eslint-disable-line @typescript-eslint/await-thenable
    await this.thread.destroy()
  }
}
