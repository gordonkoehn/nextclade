import React from 'react'

import { Progress } from 'reactstrap'
import { connect } from 'react-redux'
import i18n from 'src/i18n/i18n'

import type { State } from 'src/state/reducer'
import { AlgorithmStatus, AnylysisStatus } from 'src/state/algorithm/algorithm.state'

export function selectStatus(state: State) {
  const statusGlobal = state.algorithm.status
  const sequenceStatuses = state.algorithm.results.map(({ seqName, status }) => ({ seqName, status }))

  const parseDonePercent = 10
  let statusText = 'Idling'
  let percent = 0

  if (statusGlobal === AlgorithmStatus.parsingStarted) {
    percent = parseDonePercent
    statusText = i18n.t('Parsing...')
  } else if (statusGlobal === AlgorithmStatus.analysisStarted) {
    const total = sequenceStatuses.length
    const done = sequenceStatuses.filter(({ status }) => status === AnylysisStatus.done).length
    percent = parseDonePercent + (done / total) * (100 - parseDonePercent)
    statusText = i18n.t('Analysing sequences: {{done}}/{{total}}', { done, total })
  } else if (statusGlobal === AlgorithmStatus.done) {
    percent = 100
    statusText = i18n.t('Done')
  }

  return { percent, statusText }
}

export interface SequenceStatus {
  seqName: string
  status: AnylysisStatus
}

export interface ResultsStatusProps {
  status: { percent: number; statusText: string }
}

const mapStateToProps = (state: State) => ({
  status: selectStatus(state),
})

const mapDispatchToProps = {}

export const ResultsStatus = connect(mapStateToProps, mapDispatchToProps)(ResultsStatusDisconnected)

export function ResultsStatusDisconnected({ status }: ResultsStatusProps) {
  const { percent, statusText } = status
  const color = percent === 100 ? 'transparent' : undefined
  return (
    <div>
      <div className="text-right">{statusText}</div>
      <Progress color={color} value={percent} />
    </div>
  )
}
